A basic example on a dataset with 10 controls and 10 cases.

### 1) Count 31-mers and discard those with an abundance less than 2.


```bash
kmdiff count --file fof.txt --run-dir kcount_dir --kmer-size 31 --hard-min 2
```

### 2) Detect significant k-mers

```bash
kmdiff diff --km-run kcount_dir -1 10 -2 10 --output-dir out -s 0.01
```

Outputs:  `out/control_kmers.fasta` and `out/case_kmers.fasta`.




