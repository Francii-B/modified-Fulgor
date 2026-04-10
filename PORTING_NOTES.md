# Porting Notes

This fork starts from upstream `jermp/fulgor` tag `v4.0.0` and ports only the behavior required by the public `Phylign-Fulgor` query workflow.

## Phylign-Relevant Port

| feature | old behavior (`Francii-B/modified-Fulgor`) | new implementation point | status |
| --- | --- | --- | --- |
| threshold-union denominator | `min_score` is computed from **all** query k-mers, not only positive k-mers | [`src/ps_threshold_union.cpp`](modified-Fulgor/src/ps_threshold_union.cpp) in `pseudoalign_threshold_union_impl()` | preserved exactly |
| `--threshold 0` | treated as `min_score = 1`, so only references sharing at least one k-mer are returned | [`src/ps_threshold_union.cpp`](modified-Fulgor/src/ps_threshold_union.cpp) in `pseudoalign_threshold_union_impl()` plus CLI acceptance in [`tools/pseudoalign.cpp`](modified-Fulgor/tools/pseudoalign.cpp) | preserved exactly |
| COBS-like pseudoalign output | `pseudoalign --threshold ... --cobs` emits `*query<TAB>count` blocks followed by `_<stem><TAB><shared-kmers>` lines sorted by descending shared-k-mer count | scored threshold-union overload in [`include/index.hpp`](modified-Fulgor/include/index.hpp) and [`src/ps_threshold_union.cpp`](modified-Fulgor/src/ps_threshold_union.cpp), formatting in [`tools/pseudoalign.cpp`](modified-Fulgor/tools/pseudoalign.cpp) | preserved exactly on the Phylign-covered path |
| Phylign CLI spelling | `pseudoalign` accepts `--threshold` | argv normalization shim in [`tools/pseudoalign.cpp`](modified-Fulgor/tools/pseudoalign.cpp) | preserved exactly |
| upstream v4 output for non-COBS pseudoalign | default v4 tabular `read<TAB>count<TAB>ids...` output remains unchanged when `--cobs` is not used | [`tools/pseudoalign.cpp`](/Users/karel/git/modified-Fulgor/tools/pseudoalign.cpp) | preserved exactly |

## Intentionally Not Ported

The following modified-Fulgor features were left out on purpose because they are not required by the public `Phylign-Fulgor` workflow:

- `--best_hits`
- historical v2 index compatibility
- custom clusters / custom permutation order / `external_clusters/custom_clusters.tsv`
- index-building behavior aimed at historical Phylign datasets rather than the public query workflow


## Verification Method

Each "preserved exactly" claim is validated by regression tests comparing:

- modified-Fulgor output
- new v4.0-based port output

on identical inputs.

The following aspects are checked:
- identical reference sets per query
- identical shared k-mer counts
- identical ordering (descending counts)
- identical behavior for threshold=0

Any deviation is considered a failure.
