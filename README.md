# kmdiff

![License](http://img.shields.io/:license-affero-blue.svg)
![linux](https://github.com/tlemane/kmdiff/actions/workflows/linux.yml/badge.svg)
![osx](https://github.com/tlemane/kmdiff/actions/workflows/osx.yml/badge.svg)

This is a work in progress, it is not yet fully usable.

## Dependencies

### Build dependencies

For conveniance, all kmdiff build dependencies are included in [thirdparty directory](./thirdparty/).

### Runtime depencencies

* [kmtricks](https://github.com/tlemane/kmtricks) >= 0.0.6
* [VISOR](https://davidebolo1993.github.io/visordoc/) - (only for **kmdiff popsim**)
* [bbmap](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/) - (only for **kmdiff call**)

## Installation

### 1. Conda

A complete kmdiff runtime environment can be built using conda.

```bash
git clone https://github.com/tlemane/kmdiff
cd kmdiff
conda env create -p /kmdiff-env environment.yml
conda activate ./kmdiff-env
```

If VISOR is needed (**kmdiff popsim**):
```bash
git clone https://github.com/davidebolo1993/VISOR.git
cd VISOR
pip install -r requirements.txt
python setup.py install
chmod +x ./scripts/randomregion.r
cp ./scripts/randomregion.r ../kmdiff-env/bin
```

### 2. Build from source

#### Prerequisites
* cmake >= 3.13
* gcc >= 8.1 or XCode >= 11.0 or clang >= 7
* zlib

```bash
git clone --recursive https://github.com/tlemane/kmdiff
mkdir build ; cd build
cmake .. <args>
```

* **\<args\>**
  * `-DBUILD_KMTRICKS_FULL=[ON|OFF]` : Build kmtricks binaries.
  * `-DCOMPILE_TESTS=[ON|OFF]` : Build kmdiff tests (then run tests with `ctest --verbose`)
  * `-DSTATIC_BUILD=[ON|OFF]` : Build statically linked binary.

#### Tests

### 3. Static binary

Statically linked binaries will be provided [here](https://github.com/tlemane/kmdiff/releases) when the first stable release will be available.

## Usage

**Inputs:**
  * Control read sets.
  * Case read sets.

**Outputs:**
* Over-represented k-mers in control population.
* Over-represented k-mers in case population.

```
USAGE
  kmdiff [infos|popsim|count|diff|call]

COMMANDS
  infos  - Show build infos.
  popsim - Simulate population.
  count  - Count k-mers with kmtricks.
  diff   - Differential k-mers analysis.
  call   - SVs calling from k-mers.
```

### kmdiff count - count k-mers with kmtricks

```bash
kmdiff count --file ./inputs.fof --run-dir ./km_run --threads 8
```

**inputs.fof** (controls must appear first)
```
control1: /path/to/control1_read1.fastq ; /path/to/control1_read2.fastq
control2: /path/to/control2_read1.fastq ; /path/to/control2_read2.fastq
case1: /path/to/case1_read1.fastq ; /path/to/case1_read2.fastq
case2: /path/to/case2_read1.fastq ; /path/to/case2_read2.fastq
```

```
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
    -d --run-dir             - runtime directory 
    -k --kmer-size           - size of k-mers [16, 32] {31}
    -c --count-abundance-min - min abundance for solid k-mers {1}
    -m --max-memory          - max memory per core (in mb) {4000}

  [advanced performance tweaks]
     --minimizer-type   - minimizer type (0=lexi, 1=freq) {0}
     --minimizer-size   - size of minimizer {10}
     --repartition-type - minimizer repartition (0=unordered, 1=ordered) {0}
     --nb-partitions    - number of partitions (0=auto) {4}

  [common]
    -t --threads - Number of threads. {1}
    -h --help    - Show this message and exit. [⚑]
       --version - Show version and exit. [⚑]
    -v --verbose - Verbosity level [DEBUG|INFO|WARNING|ERROR]. {WARNING}
```

For details on `kmdiff count` parameters, see [kmtricks documentation](https://github.com/tlemane/kmtricks).

### kmdiff diff - aggregate k-mers and dump the significant ones

```bash
kmdiff diff --km-run ./km_run --output-dir ./output_diff --nb-controls 2 --nb-cases 2 --threads 8
```

```
USAGE
  kmdiff diff --km-run <DIR> --nb-controls <INT> --nb-cases <INT> [-o/--output-dir <DIR>] 
              [--significance <FLOAT>] [-c/--correction <STR>] [-t/--threads <INT>] 
              [-v/--verbose <STR>] [--kff-output] [--in-memory] [-h/--help] [--version] 

OPTIONS
  [global]
       --km-run       - kmtricks run directory. 
    -o --output-dir   - output directory. {./output}
       --nb-controls  - Number of controls. 
       --nb-cases     - Number of cases. 
       --significance - Significance threshold. {0.05}
    -c --correction   - Significance correction. {bonf}
       --kff-output   - Output significant k-mers in kff format. [⚑]
       --in-memory    - Perform correction in memory. [⚑]

  [common]
    -t --threads - Number of threads. {1}
    -h --help    - Show this message and exit. [⚑]
       --version - Show version and exit. [⚑]
    -v --verbose - Verbosity level [DEBUG|INFO|WARNING|ERROR]. {WARNING}
```

### kmdiff call - map significant k-mers with bbmap and ... (wip)

```bash
kmdiff call --reference ./ref.fa --diff-dir ./output_diff
```

```
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
    -t --threads - Number of threads. {1}
    -h --help    - Show this message and exit. [⚑]
       --version - Show version and exit. [⚑]
    -v --verbose - Verbosity level [DEBUG|INFO|WARNING|ERROR]. {WARNING}
```

## Simulation

### kmdiff popsim

wip.

## Reporting an issue

If you encounter a problem, please open an issue with the return of `kmdiff infos`, as well as the content of `kmdiff-backtrace.log` if it exists.

## Contact

Téo Lemane: teo.lemane@inria.fr

Rayan Chikhi: rayan.chikhi@pasteur.fr

Pierre Peterlongo: pierre.peterlongo@inria.fr