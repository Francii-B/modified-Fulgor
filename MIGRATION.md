# Migration

## Status

The public `Phylign-Fulgor` query workflow can run unchanged at the command level against this fork:

```bash
fulgor pseudoalign --threshold <t> -i <index>.mfur -q <query.fa> --cobs -o <output>
cat <output> | ./scripts/postprocess_cobs.py -n <n>
```

`postprocess_cobs.py` does not need to change.

## What "Unchanged" Means

If this fork is dropped in where `Phylign-Fulgor` currently expects the `modified-Fulgor` binary, the public workflow command line does not need patching.

If you keep this fork in a different checkout, the only required change is the binary path in the `Snakefile`.

## Minimal Patch When Using a Different Checkout

Replace:

```bash
./external/modified-Fulgor/build/fulgor
```

with the path to this fork's built binary, for example:

```bash
/path/to/this/fork/build-port/fulgor
```

No flag or post-processing changes are required.

In `Phylign-Fulgor/Snakefile`, the minimal edit is therefore:

```diff
- ./external/modified-Fulgor/build/fulgor pseudoalign \
+ /path/to/this/fork/build-port/fulgor pseudoalign \
```
