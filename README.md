# kmdiff

![License](http://img.shields.io/:license-affero-blue.svg)
[![kmdiff](https://github.com/tlemane/kmdiff/actions/workflows/kmdiff.yml/badge.svg)](https://github.com/tlemane/kmdiff/actions/workflows/kmdiff.yml)

## Rationale

kmdiff provides differential k-mers analysis between two populations (control and case). Each population is represented by a set of short-read sequencing. Outputs are differentially represented k-mers between controls and cases.

Statistical methods used in kmdiff are from:
* Rahman, Atif, Ingileif Hallgrímsdóttir, Michael Eisen, and Lior Pachter. "Association mapping from sequencing reads using k-mers." Elife 7 (2018): e32920. https://doi.org/10.7554/eLife.32920.001
* Mehrab Z, Mobin J, Tahmid IA, Rahman A (2021) Efficient association mapping from k-mers—An application in finding sex-specific sequences. PLOS ONE 16(1): e0245058. https://doi.org/10.1371/journal.pone.0245058.
* Patterson N, Price AL, Reich D (2006) Population Structure and Eigenanalysis. PLOS Genetics 2(12): e190. https://doi.org/10.1371/journal.pgen.0020190
* Price, A., Patterson, N., Plenge, R. et al. Principal components analysis corrects for stratification in genome-wide association studies. Nat Genet 38, 904–909 (2006). https://doi.org/10.1038/ng1847

## Dependencies

### Build dependencies

* [zlib](https://zlib.net)
* [bzip2](https://www.sourceware.org/bzip2/)

 Population stratification correction needs:
  * [GSL](https://www.gnu.org/software/gsl/)
  * [OpenBLAS](https://www.openblas.net)
  * [Lapacke](https://www.netlib.org/lapack/lapacke.html)

<details><summary><strong>Ubuntu / Debian</strong></summary>

<code>
sudo apt-get install libgsl-dev libopenblas-dev liblapacke-dev libbz2-dev zlib1g-dev zlib1g
</code>

</details>

<details><summary><strong>Fedora</strong></summary>

<code>
sudo dnf install openblas openblas-devel lapack lapack-devel gsl gsl-devel bzip2-devel
</code>

</details>

<details><summary><strong>Arch</strong></summary>

<code>
sudo pacman -S lapack lapacke openblas gsl bzip2 zlib
</code>

</details>

<details><summary><strong>macOS</strong></summary>

<code>
brew install gsl lapack openblas bzip2 zlib
</code>

</details>

For convenience, all kmdiff other build dependencies are included in [thirdparty directory](./thirdparty/).

## Installation

### 1. Conda

```bash
conda create -p /kmdiff-env
conda activate ./kmdiff-env
conda install -c bioconda -c tlemane kmdiff
```

### 2. Build from source

#### Prerequisites
* cmake >= 3.13
* gcc >= 8.1 or XCode >= 11.0 or clang >= 7
* zlib
* bzip2
* GSL + Lapacke + OpenBLAS (only with `-p 1`, [see build script](#Build-script)).

#### Clone

`git clone --recursive https://github.com/tlemane/kmdiff`

#### Build script

```
kmdiff build script.
Usage:
  ./install.sh [-r str] [-k LIST[int]] [-t int] [-c int] [-j int] [-s int] [-p] [-e] [-d] [-h]
Options:
  -r <Release|Debug> -> build type {Release}.
  -k <LIST[INT]>     -> k-mer size {"32 64 96 128"}.
  -t <0|1|2>         -> tests: 0 = disabled, 1 = compile, 2 = compile and run {2}.
  -c <1|2|4>         -> byte per count {4}.
  -j <INT>           -> nb threads {8}.
  -s <0|1>           -> population stratification correction 0 = disabled, 1 = enabled {1}
                        (-s 1 requires GSL + lapacke + OpenBLAS)
  -p                 -> compile with plugins support {disabled}
  -e                 -> use conda to install compilers/dependencies {disabled}
  -d                 -> delete cmake cache {disabled}
  -h                 -> show help.
```

If you are unable to install the prerequisites on your system, use `-e`. Compilers and build dependencies will thus be provided by a conda environment.

## Usage

### 1) `kmdiff count` - count k-mers with kmtricks

**Input file: one sample per line** (controls must appear first)
```
control1: /path/to/control1_read1.fastq ; /path/to/control1_read2.fastq
control2: /path/to/control2_read1.fastq ; /path/to/control2_read2.fastq
case1: /path/to/case1_read1.fastq ; /path/to/case1_read2.fastq
case2: /path/to/case2_read1.fastq ; /path/to/case2_read2.fastq
```

Supported files: fasta/fastq, gzipped or not.

**Options**
```
kmdiff count v1.0.0

DESCRIPTION
  Count k-mers with kmtricks.

USAGE
  kmdiff count -f/--file <FILE> -d/--run-dir <DIR> [-k/--kmer-size <INT>] [-c/--hard-min <INT>]
               [-r/--recurrence-min <INT>] [--minimizer-type <INT>]
               [--minimizer-size <INT>] [--repartition-type <INT>] [--nb-partitions <INT>]
               [-t/--threads <INT>] [-v/--verbose <STR>] [-h/--help] [--version]

OPTIONS
  [global]
    -f --file           - fof that contains path of read files
    -d --run-dir        - output directory.
    -k --kmer-size      - size of k-mers [8, 127] {31}
    -c --hard-min       - min abundance to keep a k-mer {1}
    -r --recurrence-min - min recurrence to keep a k-mer {1}

  [advanced performance tweaks]
     --minimizer-type   - minimizer type (0=lexi, 1=freq) {0}
     --minimizer-size   - size of minimizer {10}
     --repartition-type - minimizer repartition (0=unordered, 1=ordered) {0}
     --nb-partitions    - number of partitions (0=auto) {0}

  [common]
    -t --threads - number of threads. {8}
    -h --help    - show this message and exit. [⚑]
       --version - show version and exit. [⚑]
    -v --verbose - Verbosity level [debug|info|warning|error]. {info}
```

### 2) `kmdiff diff` - aggregate k-mers and dump the significant ones

```
kmdiff diff v1.0.0

DESCRIPTION
  Differential k-mers analysis.

USAGE
  kmdiff diff -d/--km-run <DIR> -1/--nb-controls <INT> -2/--nb-cases <INT> [-o/--output-dir <DIR>]
              [-s/--significance <FLOAT>] [-u/--cutoff <INT>] [-c/--correction <STR>]
              [--gender <FILE>] [--kmer-pca <FLOAT>] [--ploidy <INT>] [--n-pc <INT>]
              [-t/--threads <INT>] [-v/--verbose <STR>] [-f/--kff-output] [-m/--in-memory]
              [--keep-tmp] [--pop-correction] [-h/--help] [--version]

OPTIONS
  [global]
    -d --km-run       - kmtricks run directory.
    -o --output-dir   - output directory. {./kmdiff_output}
    -1 --nb-controls  - number of controls.
    -2 --nb-cases     - number of cases.
    -s --significance - significance threshold. {0.05}
    -u --cutoff       - Divide the significance threshold by N.
                        Since a large number of k-mers are tested, k-mers with p-values too close to the significance
                        threshold will not pass the last steps of correction.
                        It allows to discard some k-mers a bit earlier and thus save space and time. {100000}
    -c --correction   - significance correction. (bonferroni|benjamini|sidak|holm|disabled) {bonferroni}
    -f --kff-output   - output significant k-mers in kff format. [⚑]
    -m --in-memory    - in-memory correction. [⚑]
       --keep-tmp     - keep tmp files. [⚑]

  [population stratification]
     --pop-correction - apply correction for population stratification. [⚑]
     --gender         - gender file, one sample per line with the id and the gender (M,F,U), space-separated.
     --kmer-pca       - proportion of k-mers used for PCA (in [0.0, 0.05]). {0.001}
     --ploidy         - ploidy level. {2}
     --n-pc           - number of principal components (in [2, 10]). {2}

  [common]
    -t --threads - number of threads. {8}
    -h --help    - show this message and exit. [⚑]
       --version - show version and exit. [⚑]
    -v --verbose - Verbosity level [debug|info|warning|error]. {info}
```

**Outputs**
* control significant k-mers: `<output_dir>/control_kmers.[fasta|kff]`
* case significant k-mers: `<output_dir>/case_kmers.[fasta|kff]`

Abundances and p-values are provided in fasta headers.

## Testing

An example on a small dataset is available [here](./examples).

## Reporting an issue

If you encounter a problem, please open an issue with the return of `kmdiff infos`, as well as the content of `kmdiff-backtrace.log` if it exists.

## Contact

Téo Lemane: teo[dot]lemane[at]inria[dot]fr  
Rayan Chikhi: rayan[dot]chikhi[at]pasteur[dot]fr  
Pierre Peterlongo: pierre[dot]peterlongo[at]inria[dot]fr  
