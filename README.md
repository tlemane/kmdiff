# kmdiff

![License](http://img.shields.io/:license-affero-blue.svg)
![linux](https://github.com/tlemane/kmdiff/actions/workflows/linux.yml/badge.svg)
![osx](https://github.com/tlemane/kmdiff/actions/workflows/osx.yml/badge.svg)

This is a work in progress, it is not yet fully usable.

## Rationale

kmdiff takes two read sets and outputs differentially represented k-mers between controls and cases.

Statistical methods used in kmdiff are from:
* Rahman, Atif, Ingileif Hallgrímsdóttir, Michael Eisen, and Lior Pachter. "Association mapping from sequencing reads using k-mers." Elife 7 (2018): e32920. https://doi.org/10.7554/eLife.32920.001
* Mehrab Z, Mobin J, Tahmid IA, Rahman A (2021) Efficient association mapping from k-mers—An application in finding sex-specific sequences. PLOS ONE 16(1): e0245058. https://doi.org/10.1371/journal.pone.0245058.
* Patterson N, Price AL, Reich D (2006) Population Structure and Eigenanalysis. PLOS Genetics 2(12): e190. https://doi.org/10.1371/journal.pgen.0020190
* Price, A., Patterson, N., Plenge, R. et al. Principal components analysis corrects for stratification in genome-wide association studies. Nat Genet 38, 904–909 (2006). https://doi.org/10.1038/ng1847

## Dependencies

### Build dependencies


A thirdparty tool used for population stratification correction needs:
  * [GSL](https://www.gnu.org/software/gsl/)
  * [OpenBLAS](https://www.openblas.net)
  * [Lapack](http://www.netlib.org/lapack/)

For conveniance, all kmdiff other build dependencies are included in [thirdparty directory](./thirdparty/).

### Runtime dependencies

Runtime dependencies are included in the conda recipes.

If you don't use conda, the following tools must be in your `PATH` (only for `kmdiff call`):
  * [samtools](http://www.htslib.org)
  * [BBMap](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/)
  * [ABySS](https://github.com/bcgsc/abyss)

## Installation

### 1. Conda

```bash
conda create -p /kmdiff-env
conda activate ./kmdiff-env
conda install -c bioconda -c tlemane -c tlemane/label/dev kmdiff # temporary
```

### 2. Build from source

#### Prerequisites
* cmake >= 3.13
* gcc >= 8.1 or XCode >= 11.0 or clang >= 7
* zlib

#### Clone

`git clone --recursive https://github.com/tlemane/kmdiff`

#### Build script

```
kmdiff build script.
Usage: 
  ./install.sh [-r str] [-k int] [-t int] [-c int] [-p int] [-j int] [-d] [-s] [-h]
Options: 
  -r <Release|Debug> -> build type {Release}.
  -k <0|32|64>       -> max k-mer size, 0 = both {0}.
  -t <0|1|2>         -> tests: 0 = disable, 1 = compile, 2 = compile and run {2}.
  -c <1|2|4>         -> byte per count {4}.
  -j <INT>           -> nb threads {8}.
  -p <0|1>           -> population stratification correction 0 = disable, 1 = enable {1}
                        (-p 1 requires GSL + lapack + OpenBLAS)
  -d                 -> dev build {disable}.
  -s                 -> static build {disable}.
  -h                 -> show help.
```

`-k 0` produces two binaries:
* `kmdiff` supports k in [16, 32]
* `kmdiff64` supports k in [32, 64]

### 3. Static binary

Statically linked binaries will be provided [here](https://github.com/tlemane/kmdiff/releases) when the first stable release will be available.

## Usage

```
USAGE
  kmdiff [infos|count|diff|call]

COMMANDS
  infos  - Show build infos.
  count  - Count k-mers with kmtricks.
  diff   - Differential k-mers analysis.
  call   - SVs calling from k-mers (WIP).
```

### 1) kmdiff count - count k-mers with kmtricks

**Input file: one sample per line** (controls must appear first)
```
control1: /path/to/control1_read1.fastq ; /path/to/control1_read2.fastq
control2: /path/to/control2_read1.fastq ; /path/to/control2_read2.fastq
case1: /path/to/case1_read1.fastq ; /path/to/case1_read2.fastq
case2: /path/to/case2_read1.fastq ; /path/to/case2_read2.fastq
```

**A sample-specific min abundance threshold can be specified as follows:**
```
control1: /path/to/control1_read1.fastq ; /path/to/control1_read2.fastq ! 3
control2: /path/to/control2_read1.fastq ; /path/to/control2_read2.fastq ! 5
case1: /path/to/case1_read1.fastq ; /path/to/case1_read2.fastq ! 2
case2: /path/to/case2_read1.fastq ; /path/to/case2_read2.fastq ! 2
```

Supported files: fasta/fastq, gzipped or not.

**Usage**
```bash
kmdiff count --file ./inputs.fof --run-dir ./km_run --threads 8
```

**Options**
```
kmdiff count v0.1

DESCRIPTION
  Count k-mers with kmtricks.

USAGE
  kmdiff count -f/--file <FILE> -d/--run-dir <DIR> [-k/--kmer-size <INT>] [-c/--count-abundance-min <INT>] 
               [-m/--max-memory <INT>] [--minimizer-type <INT>] [--minimizer-size <INT>] 
               [--repartition-type <INT>] [--nb-partitions <INT>] [-t/--threads <INT>] 
               [-v/--verbose <STR>] [-h/--help] [--version] 

OPTIONS
  [global]
    -f --file                - fof that contains path of read files 
    -d --run-dir             - Output directory. 
    -k --kmer-size           - size of k-mers [16, 32] {31}
    -c --count-abundance-min - min abundance for solid k-mers {1}
    -m --max-memory          - max memory per core (in mb) {4000}

  [advanced performance tweaks]
     --minimizer-type   - minimizer type (0=lexi, 1=freq) {0}
     --minimizer-size   - size of minimizer {10}
     --repartition-type - minimizer repartition (0=unordered, 1=ordered) {0}
     --nb-partitions    - number of partitions (0=auto) {4}

  [common]
    -t --threads - Number of threads. {8}
    -h --help    - Show this message and exit. [⚑]
       --version - Show version and exit. [⚑]
    -v --verbose - Verbosity level [DEBUG|INFO|WARNING|ERROR]. {INFO}
```

### 2) kmdiff diff - aggregate k-mers and dump the significant ones

**Usage**
```bash
kmdiff diff --km-run ./km_run --output-dir ./output_diff --nb-controls 2 --nb-cases 2 --threads 8
```

**Options**
```
kmdiff diff v0.1

DESCRIPTION
  Differential k-mers analysis.

USAGE
  kmdiff diff --km-run <DIR> --nb-controls <INT> --nb-cases <INT> [-o/--output-dir <DIR>] 
              [--coverage <INT>] [--significance <FLOAT>] [-c/--correction <STR>] 
              [--kmer-pca <FLOAT>] [--ploidy <INT>] [--n-pc <INT>] [-t/--threads <INT>] 
              [-v/--verbose <STR>] [--kff-output] [--in-memory] [--pop-correction] 
              [-h/--help] [--version] 

OPTIONS
  [global]
       --km-run       - kmtricks run directory. 
    -o --output-dir   - output directory. {./kmdiff_output}
       --nb-controls  - Number of controls. 
       --nb-cases     - Number of cases. 
       --coverage     - Coverage (running time concern, no impact on results). {20}
       --significance - Significance threshold. {0.05}
    -c --correction   - Significance correction. {bonferroni}
       --kff-output   - Output significant k-mers in kff format. [⚑]
       --in-memory    - Perform correction in memory. [⚑]

  [population stratification] - (WIP)
     --pop-correction - Apply correction of population stratification. [⚑]
     --kmer-pca       - Proportion of k-mers used for PCA. {0.001}
     --ploidy         - Ploidy level. {2}
     --n-pc           - Number of principal components. {2}

  [common]
    -t --threads - Number of threads. {8}
    -h --help    - Show this message and exit. [⚑]
       --version - Show version and exit. [⚑]
    -v --verbose - Verbosity level [DEBUG|INFO|WARNING|ERROR]. {INFO}
```

**Outputs**
* control significant k-mers: `<output_dir>/control_kmers.[fasta|kff]`
* case significant k-mers: `<output_dir>/case_kmers.fasta.[fasta|kff]`

**Warning**: The next module `kmdiff call` doesn't support kff format.

### 3) kmdiff call - (WIP)

```bash
kmdiff call --reference ./ref.fa --diff-dir ./output_diff
```

```
kmdiff call v0.1

DESCRIPTION
  SVs calling from k-mers.

USAGE
  kmdiff call -r/--reference <FILE> -d/--diff-dir <DIR> [--seed-size <INT>] [-t/--threads <INT>] 
              [-v/--verbose <STR>] [-h/--help] [--version] 

OPTIONS
  [global]
    -r --reference - Reference genome. 
    -d --diff-dir  - Output directory of kmdiff diff. 
       --seed-size - Size of seeds for mapping. {10}

  [common]
    -t --threads - Number of threads. {8}
    -h --help    - Show this message and exit. [⚑]
       --version - Show version and exit. [⚑]
    -v --verbose - Verbosity level [DEBUG|INFO|WARNING|ERROR]. {INFO}
```

## Reporting an issue

If you encounter a problem, please open an issue with the return of `kmdiff infos`, as well as the content of `kmdiff-backtrace.log` if it exists.

## Contact

Téo Lemane: teo[dot]lemane[at]inria[dot]fr  
Rayan Chikhi: rayan[dot]chikhi[at]pasteur[dot]fr  
Pierre Peterlongo: pierre[dot]peterlongo[at]inria[dot]fr  