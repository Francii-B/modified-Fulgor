# Porting Notes

This fork starts from upstream `jermp/fulgor` tag `v4.0.0` and ports only the behavior required by the public `Phylign-Fulgor` query workflow.

## Design Constraint

Only behavior required by the public `Phylign-Fulgor` query workflow is ported here. Unrelated historical or custom `modified-Fulgor` features were intentionally left out so the fork stays small, reviewable, and aligned with upstream `v4.0.0`.

## Phylign-Relevant Port

| feature | old behavior (`Francii-B/modified-Fulgor`) | new implementation point | status |
| --- | --- | --- | --- |
| threshold-union denominator | `min_score` is computed from **all** query k-mers, not only positive k-mers | `src/ps_threshold_union.cpp` in `pseudoalign_threshold_union_impl()` | preserved exactly for the covered query path |
| `--threshold 0` | treated as `min_score = 1`, so only references sharing at least one k-mer are returned | `src/ps_threshold_union.cpp` in `pseudoalign_threshold_union_impl()` plus CLI acceptance in `tools/pseudoalign.cpp` | preserved exactly for the covered query path |
| COBS-like pseudoalign output | `pseudoalign --threshold ... --cobs` emits `*query<TAB>count` blocks followed by `_<stem><TAB><shared-kmers>` lines sorted by descending shared-k-mer count | scored threshold-union overload in `include/index.hpp` and `src/ps_threshold_union.cpp`, formatting in `tools/pseudoalign.cpp` | preserved exactly for the covered Phylign path |
| Phylign CLI spelling | `pseudoalign` accepts `--threshold` | argv normalization shim in `tools/pseudoalign.cpp` | preserved exactly for the covered public command path |
| upstream v4 output for non-COBS pseudoalign | default v4 tabular `read<TAB>count<TAB>ids...` output remains unchanged when `--cobs` is not used | `tools/pseudoalign.cpp` | intentionally kept on the upstream v4 code path |

## Differences from Upstream v4.0.0

Relative to vanilla upstream `jermp/fulgor` tag `v4.0.0`, this fork changes only the Phylign-relevant `pseudoalign` surface:

- threshold denominator semantics use all query k-mers when computing the threshold-union minimum score; upstream v4 uses its own default threshold-union semantics.
- `--threshold 0` is accepted and treated as a positive-only query path (`min_score = 1`); upstream v4 rejects `0` for the threshold option.
- `--cobs` is available on `pseudoalign` and emits the Phylign-compatible block format consumed by `scripts/postprocess_cobs.py`; upstream v4 does not provide that Phylign-oriented output mode.
- `--threshold` is accepted as a compatibility spelling for the public Phylign command line; upstream v4 exposes the threshold option through its native CLI spelling.

## Intentionally Not Ported

The following modified-Fulgor features were left out on purpose because they are not required by the public `Phylign-Fulgor` workflow:

- `--best_hits`
- historical v2 index compatibility
- custom clusters / custom permutation order / `external_clusters/custom_clusters.tsv`
- index-building behavior aimed at historical Phylign datasets rather than the public query workflow


## Verification Method

Each "preserved exactly" claim is validated by regression tests comparing:

- the old `modified-Fulgor` behavior
- this v4.0-based port

on identical inputs.

The following aspects are checked:
- identical reference sets per query
- identical shared-k-mer counts
- identical ordering when ordering is part of the legacy contract
- identical `--threshold 0` behavior

Any deviation is considered a failure.

## Test Coverage

`tests/phylign_port_regression.py` is the regression contract for this port. It is expected to validate:

- `pseudoalign` on the Phylign command path with `--cobs`
- threshold behavior across representative thresholds, including a non-zero threshold case that distinguishes upstream v4 from the legacy behavior
- `--threshold 0` behavior on the legacy-compatible path
- compatibility of `--cobs` output with `scripts/postprocess_cobs.py`
- cross-checks against the legacy `modified-Fulgor` binary for identical reference sets, shared-k-mer counts, and ordering where ordering is part of the legacy contract
