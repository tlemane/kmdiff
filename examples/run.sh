#!/bin/bash

kmdiff count --file fof.txt --run-dir kc_dir --kmer-size 31 --hard-min 2

kmdiff diff --km-run kc_dir -1 10 -2 10 --output-dir out -s 0.01
