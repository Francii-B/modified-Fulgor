#!/usr/bin/env python3

import argparse
import os
import subprocess
import sys
import tempfile
from pathlib import Path


EXPECTED = {
    "threshold_normal_upstream": "q_threshold\t2\t0\t2\n",
    "threshold_normal_old_new": "q_threshold\t0\n",
    "threshold_zero_old_new": "q_threshold\t3\t0\t1\t2\n",
    "cobs_old_new": "*q_cobs\t3\n_001_alpha\t3\n_003_gamma\t3\n_002_beta\t2\n",
    "postprocess_old_new": "*q_cobs\t3\n_001_alpha\t3\n_003_gamma\t3\n",
}


def run(cmd, *, input_text=None, expect_success=True):
    result = subprocess.run(
        cmd,
        input=input_text,
        text=True,
        capture_output=True,
        check=False,
    )
    if expect_success and result.returncode != 0:
        raise AssertionError(
            f"command failed ({result.returncode}): {' '.join(cmd)}\n"
            f"stdout:\n{result.stdout}\n"
            f"stderr:\n{result.stderr}"
        )
    return result


def build_indexes(binary, refs_dir, work_dir):
    refs = sorted(refs_dir.glob("*.fa"))
    list_path = work_dir / "filenames.txt"
    list_path.write_text("".join(f"{path.resolve()}\n" for path in refs), encoding="utf-8")

    tmp_dir = work_dir / "tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    out_base = work_dir / "tiny"
    cmd = [
        str(binary),
        "build",
        "-l",
        str(list_path),
        "-o",
        str(out_base),
        "-k",
        "5",
        "-m",
        "3",
        "-d",
        str(tmp_dir),
        "-g",
        "1",
        "-t",
        "1",
        "--meta",
    ]
    run(cmd)
    fur = out_base.with_suffix(".fur")
    mfur = out_base.with_suffix(".mfur")
    if not fur.exists() or not mfur.exists():
        raise AssertionError(f"expected indexes were not created under {work_dir}")
    return fur, mfur


def run_pseudoalign(binary, index_path, query_path, output_path, *, threshold_flag, threshold,
                    cobs=False, expect_success=True):
    cmd = [
        str(binary),
        "pseudoalign",
        threshold_flag,
        threshold,
        "-t",
        "1",
        "-i",
        str(index_path),
        "-q",
        str(query_path),
    ]
    if cobs:
        cmd.append("--cobs")
    cmd.extend(["-o", str(output_path)])
    return run(cmd, expect_success=expect_success)


def assert_file_text(path, expected):
    got = path.read_text(encoding="utf-8")
    if got != expected:
        raise AssertionError(f"unexpected contents for {path}\nexpected:\n{expected}\ngot:\n{got}")


def run_postprocess(script_path, cobs_output):
    cmd = [sys.executable, str(script_path), "-n", "1"]
    return run(cmd, input_text=cobs_output)


def exercise_binary(binary, label, refs_dir, query_threshold, query_cobs, postprocess_script):
    with tempfile.TemporaryDirectory(prefix=f"phylign-port-{label}-") as tmp:
        work_dir = Path(tmp)
        fur, mfur = build_indexes(binary, refs_dir, work_dir)

        threshold_out = work_dir / "threshold.tsv"
        threshold_flag = "-r" if label == "upstream" else "--threshold"
        run_pseudoalign(binary, fur, query_threshold, threshold_out, threshold_flag=threshold_flag,
                        threshold="0.75")
        expected_threshold = (
            EXPECTED["threshold_normal_upstream"]
            if label == "upstream"
            else EXPECTED["threshold_normal_old_new"]
        )
        assert_file_text(threshold_out, expected_threshold)

        zero_out = work_dir / "threshold_zero.tsv"
        if label == "upstream":
            result = run_pseudoalign(
                binary,
                fur,
                query_threshold,
                zero_out,
                threshold_flag="-r",
                threshold="0",
                expect_success=False,
            )
            if result.returncode == 0:
                raise AssertionError("upstream v4 unexpectedly accepted -r 0")
        else:
            run_pseudoalign(binary, fur, query_threshold, zero_out, threshold_flag="--threshold",
                            threshold="0")
            assert_file_text(zero_out, EXPECTED["threshold_zero_old_new"])

        cobs_out = work_dir / "cobs.tsv"
        if label == "upstream":
            result = run_pseudoalign(
                binary,
                mfur,
                query_cobs,
                cobs_out,
                threshold_flag="--threshold",
                threshold="0",
                cobs=True,
                expect_success=False,
            )
            if result.returncode == 0:
                raise AssertionError("upstream v4 unexpectedly accepted the Phylign COBS command")
            return

        run_pseudoalign(binary, mfur, query_cobs, cobs_out, threshold_flag="--threshold",
                        threshold="0", cobs=True)
        assert_file_text(cobs_out, EXPECTED["cobs_old_new"])

        if postprocess_script is not None:
            postprocess = run_postprocess(postprocess_script, cobs_out.read_text(encoding="utf-8"))
            if postprocess.stdout != EXPECTED["postprocess_old_new"]:
                raise AssertionError(
                    "unexpected postprocess_cobs.py output\n"
                    f"expected:\n{EXPECTED['postprocess_old_new']}\n"
                    f"got:\n{postprocess.stdout}"
                )


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Build tiny .fur/.mfur fixtures and compare upstream v4.0.0, old modified-Fulgor, "
            "and the ported binary on the Phylign-relevant pseudoalign surface."
        )
    )
    parser.add_argument("--new-binary", required=True, type=Path)
    parser.add_argument("--old-binary", type=Path, default=os.environ.get("FULGOR_OLD_BINARY"))
    parser.add_argument(
        "--upstream-binary", type=Path, default=os.environ.get("FULGOR_UPSTREAM_BINARY")
    )
    parser.add_argument(
        "--postprocess-script",
        type=Path,
        default=os.environ.get("PHYALIGN_POSTPROCESS_COBS"),
    )
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parent
    data_dir = repo_root / "data" / "phylign_port"
    refs_dir = data_dir / "refs"
    query_threshold = data_dir / "query_threshold.fa"
    query_cobs = data_dir / "query_cobs.fa"

    exercise_binary(args.new_binary, "new", refs_dir, query_threshold, query_cobs,
                    args.postprocess_script)

    if args.old_binary is not None:
        exercise_binary(Path(args.old_binary), "old", refs_dir, query_threshold, query_cobs,
                        args.postprocess_script)

    if args.upstream_binary is not None:
        exercise_binary(Path(args.upstream_binary), "upstream", refs_dir, query_threshold,
                        query_cobs, args.postprocess_script)

    if args.old_binary is None or args.upstream_binary is None:
        print(
            "new-binary regression checks passed; full cross-version comparison skipped because "
            "--old-binary and/or --upstream-binary were not provided",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
