#!/usr/bin/env bash

set -euo pipefail

usage() {
  cat <<'EOF'
Publish conda packages to the gh-pages branch as a static conda channel.

Usage:
  publish_conda_channel.sh [options] PKG [PKG ...]
  publish_conda_channel.sh [options] --from-conda-bld DIR

Options:
  -b, --branch BRANCH        Channel branch to publish to [default: gh-pages]
  -r, --remote REMOTE        Git remote to push to [default: origin]
  -m, --message MESSAGE      Commit message
      --subdir SUBDIR        Force a subdir for explicitly listed packages
      --from-conda-bld DIR   Import all *.tar.bz2 and *.conda from DIR/<subdir>/
      --no-push              Commit locally but do not push
      --keep-worktree        Do not delete the temporary worktree
  -h, --help                 Show this help

Examples:
  publish_conda_channel.sh ~/miniconda/conda-bld/osx-arm64/modified-fulgor-2.0.0.post2-hcb8d3e5_1.tar.bz2

  publish_conda_channel.sh --from-conda-bld ~/miniconda/conda-bld

  publish_conda_channel.sh --subdir linux-64 /tmp/modified-fulgor-2.0.0.post2-*.conda
EOF
}

BRANCH="gh-pages"
REMOTE="origin"
MESSAGE=""
FORCE_SUBDIR=""
FROM_CONDA_BLD=""
NO_PUSH=0
KEEP_WORKTREE=0

ARGS=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    -b|--branch)
      BRANCH="$2"
      shift 2
      ;;
    -r|--remote)
      REMOTE="$2"
      shift 2
      ;;
    -m|--message)
      MESSAGE="$2"
      shift 2
      ;;
    --subdir)
      FORCE_SUBDIR="$2"
      shift 2
      ;;
    --from-conda-bld)
      FROM_CONDA_BLD="$2"
      shift 2
      ;;
    --no-push)
      NO_PUSH=1
      shift
      ;;
    --keep-worktree)
      KEEP_WORKTREE=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    --)
      shift
      while [[ $# -gt 0 ]]; do
        ARGS+=("$1")
        shift
      done
      ;;
    -*)
      echo "Unknown option: $1" >&2
      usage >&2
      exit 2
      ;;
    *)
      ARGS+=("$1")
      shift
      ;;
  esac
done

if ! git rev-parse --show-toplevel >/dev/null 2>&1; then
  echo "Error: run this from inside the git repository." >&2
  exit 1
fi

REPO_ROOT="$(git rev-parse --show-toplevel)"
cd "$REPO_ROOT"

if [[ -n "$FROM_CONDA_BLD" && ${#ARGS[@]} -gt 0 ]]; then
  echo "Error: use either explicit package paths or --from-conda-bld, not both." >&2
  exit 1
fi

if [[ -z "$FROM_CONDA_BLD" && ${#ARGS[@]} -eq 0 ]]; then
  echo "Error: no packages specified." >&2
  usage >&2
  exit 1
fi

if ! python - <<'PY' >/dev/null 2>&1
import importlib.util, sys
sys.exit(0 if importlib.util.find_spec("conda_index") else 1)
PY
then
  echo "Error: python package 'conda_index' is not installed." >&2
  echo "Install it with: python -m pip install conda-index" >&2
  exit 1
fi

is_subdir_name() {
  local d="$1"
  [[ "$d" =~ ^(noarch|linux-.*|osx-.*|win-.*)$ ]]
}

detect_subdir_for_file() {
  local pkg="$1"
  local parent
  parent="$(basename "$(dirname "$pkg")")"

  if [[ -n "$FORCE_SUBDIR" ]]; then
    printf '%s\n' "$FORCE_SUBDIR"
    return 0
  fi

  if is_subdir_name "$parent"; then
    printf '%s\n' "$parent"
    return 0
  fi

  return 1
}

collect_from_conda_bld() {
  local root="$1"
  local -n out_files_ref="$2"
  local d

  if [[ ! -d "$root" ]]; then
    echo "Error: conda-bld directory not found: $root" >&2
    exit 1
  fi

  while IFS= read -r -d '' d; do
    while IFS= read -r -d '' f; do
      out_files_ref+=("$f")
    done < <(find "$d" -maxdepth 1 -type f \( -name '*.tar.bz2' -o -name '*.conda' \) -print0 | sort -z)
  done < <(find "$root" -mindepth 1 -maxdepth 1 -type d -print0 | sort -z)
}

PKGS=()
if [[ -n "$FROM_CONDA_BLD" ]]; then
  collect_from_conda_bld "$FROM_CONDA_BLD" PKGS
else
  PKGS=("${ARGS[@]}")
fi

if [[ ${#PKGS[@]} -eq 0 ]]; then
  echo "Error: no package artifacts found." >&2
  exit 1
fi

for pkg in "${PKGS[@]}"; do
  if [[ ! -f "$pkg" ]]; then
    echo "Error: package file not found: $pkg" >&2
    exit 1
  fi
done

TMP_WORKTREE="$(mktemp -d "${TMPDIR:-/tmp}/conda-channel.XXXXXX")"

cleanup() {
  local status=$?
  if git worktree list --porcelain | grep -Fq "worktree $TMP_WORKTREE"; then
    git worktree remove --force "$TMP_WORKTREE" >/dev/null 2>&1 || true
  fi
  if [[ $KEEP_WORKTREE -eq 1 ]]; then
    echo "Kept temporary worktree at: $TMP_WORKTREE" >&2
  else
    rm -rf "$TMP_WORKTREE" >/dev/null 2>&1 || true
  fi
  exit $status
}
trap cleanup EXIT

if git show-ref --verify --quiet "refs/heads/$BRANCH"; then
  git worktree add --force "$TMP_WORKTREE" "$BRANCH" >/dev/null
elif git ls-remote --exit-code --heads "$REMOTE" "$BRANCH" >/dev/null 2>&1; then
  git fetch "$REMOTE" "$BRANCH:$BRANCH"
  git worktree add --force "$TMP_WORKTREE" "$BRANCH" >/dev/null
else
  git worktree add --detach "$TMP_WORKTREE" >/dev/null
  (
    cd "$TMP_WORKTREE"
    git checkout --orphan "$BRANCH" >/dev/null 2>&1
    find . -mindepth 1 -maxdepth 1 ! -name .git -exec rm -rf {} +
  )
fi

cd "$TMP_WORKTREE"

mkdir -p noarch
touch .nojekyll

ADDED=()

for pkg in "${PKGS[@]}"; do
  subdir="$(detect_subdir_for_file "$pkg")" || {
    echo "Error: cannot infer platform subdir for $pkg" >&2
    echo "Put the file under a standard conda-bld subdir or pass --subdir." >&2
    exit 1
  }

  mkdir -p "$subdir"
  cp -f "$pkg" "$subdir/"
  ADDED+=("$subdir/$(basename "$pkg")")
done

python -m conda_index "$TMP_WORKTREE"

git add -A

if git diff --cached --quiet; then
  echo "No channel changes detected."
  exit 0
fi

if [[ -z "$MESSAGE" ]]; then
  MESSAGE="Update conda channel with $(printf '%s ' "${ADDED[@]}")"
fi

git commit -m "$MESSAGE"

if [[ $NO_PUSH -eq 0 ]]; then
  git push "$REMOTE" "$BRANCH"
fi

echo
echo "Published to branch: $BRANCH"
echo "Artifacts added or updated:"
printf '  - %s\n' "${ADDED[@]}"
echo
echo "Channel root now contains:"
find . -maxdepth 2 -type f | sort
