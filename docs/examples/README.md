# Examples

## Stage 1 simulation helper

`generate_stage1_sim.py` simulates phenotypes and genotypes, runs the Stage 1
pipeline, and writes `.stage1.tsv.gz`, `.cov.gz`, and metadata for Stage 2
validation.

```bash
python docs/examples/generate_stage1_sim.py \
  --out-prefix /tmp/stage1_sim \
  --n-samples 300000 \
  --n-snps 100
```
