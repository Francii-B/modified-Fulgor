Fulgor
======

Fulgor is a *(meta-) colored compacted de Bruijn graph* index for large-scale matching and color queries, powered by [SSHash](https://github.com/jermp/sshash) and [GGCAT](https://github.com/algbio/GGCAT).

The Fulgor index is described in the following papers:

- [**Fulgor: A Fast and Compact k-mer Index for Large-Scale Matching and Color Queries**](https://almob.biomedcentral.com/articles/10.1186/s13015-024-00251-9) (Algorithms for Molecular Biology, ALMOB 2024), and
<!--	- [WABI 2023 version](https://drops.dagstuhl.de/opus/volltexte/2023/18644/)-->

- [**Meta-colored compacted de Bruijn graphs**](https://www.biorxiv.org/content/10.1101/2023.07.21.550101v2.full.pdf) (International Conference on Research in Computational Molecular Biology, RECOMB 2024).

Please, cite these papers if you use Fulgor.

### Table of contents
* [Dependencies](#dependencies)
* [Compiling the code](#compiling-the-code)
* [Tools and usage](#tools-and-usage)
* [Quick start](#quick-start)
* [Indexing an example Salmonella pangenome](#indexing-an-example-salmonella-pangenome)
* [Pseudoalignment output format](#pseudoalignment-output-format)

Dependencies
------------

#### GGCAT

The code uses the [GGCAT](https://github.com/algbio/GGCAT) Rust library,
so make sure you have Rust installed. If not, Rust can be installed as recommended [here](https://www.rust-lang.org/tools/install), with

	curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

#### zlib

If you do not have `zlib` installed, you can do

    sudo apt-get install zlib1g

if you are on Linux/Ubuntu, or

    brew install zlib

if you are using MacOS.


Compiling the code
------------------

The code is tested on Linux with `gcc` and on MacOS with `clang`.
To build the code, [`CMake`](https://cmake.org/) is required.

First clone the repository with

    git clone https://github.com/jermp/fulgor.git

and then do

    git submodule update --init --recursive

to pull all necessary submodules before compilation.

To compile the code for a release environment (see file `CMakeLists.txt` for the used compilation flags), it is sufficient to do the following, within the parent `fulgor` directory:

    mkdir build
    cd build
    cmake ..
    make -j

For a testing environment, use the following instead:

    mkdir debug_build
    cd debug_build
    cmake .. -D CMAKE_BUILD_TYPE=Debug -D FULGOR_USE_SANITIZERS=On
    make -j


Tools and usage
---------------

There is one executable called `fulgor` after the compilation, which can be used to run a tool.
Run `./fulgor` to see a list of available tools.

	== Fulgor: a (meta-) colored compacted de Bruijn graph index =============================

	Usage: ./fulgor <tool> ...

	Tools:
	  build                  build a Fulgor index
	  pseudoalign            pseudoalign reads to references
	  stats                  print index statistics
	  print-filenames        print all reference filenames

	Advanced tools:
	  meta                   partition a Fulgor index and build a meta-colored Fulgor index
	  differential           partition a Fulgor index and build a differential-colored Fulgor index
	  meta-differential      partition a meta-Fulgor index and build a meta-differential-colored Fulgor index
	  dump                   write colors to an output file in text format

For large-scale indexing, it could be necessary to increase the number of file descriptors that can be opened simultaneously:

	ulimit -n 2048

Quick start
-----------

This short demo shows how to index the 10-genome collection
in the folder `test_data/salmonella_10` with Fulgor.
We will use the standard value k = 31.

First create a list of filenames (with absolute paths) for the files in `test_data/salmonella_10`.
From `fulgor/test_data`, do

	find $(pwd)/salmonella_10/* > salmonella_10_filenames.txt

Then, from `fulgor/build`, run

	./fulgor build -l ../test_data/salmonella_10_filenames.txt -o ../test_data/salmonella_10 -k 31 -m 19 -d tmp_dir -g 1 -t 1 --verbose --check

to build an index that will be serialized to the file `test_data/salmonella_10.fur`.


Indexing an example Salmonella pangenome
----------------------------------------

In this example, we will build a Fulgor index, with k = 31, for the 4,546 Salmonella genomes that can be downloaded from [here](https://zenodo.org/record/1323684).

We assume all commands are issue from within the home (`~/`) directory.

After download,
create a list of all `.fasta` filenames with

	find $(pwd)/Salmonella_enterica/Genomes/*.fasta > salmonella_4546_filenames.txt

and, from `fulgor/build`, run

	./fulgor build -l ~/salmonella_4546_filenames.txt -o ~/Salmonella_enterica/salmonella_4546 -k 31 -m 20 -d tmp_dir -g 8 -t 8 --verbose --check

which will create an index named `~/Salmonella_enterica/salmonella_4546.fur` of 0.266 GB.

We can now pseudoalign the reads from SRR801268, as follows.

First, download the reads in `~/` with (assuming you have `wget` installed):

	cd
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR801/SRR801268/SRR801268_1.fastq.gz

and then process them with:

	./fulgor pseudoalign -i ~/Salmonella_enterica/salmonella_4546.fur -q ~/SRR801268_1.fastq.gz -t 8 -o /dev/null

	mapped 6584304 reads
	elapsed = 130133 millisec / 130.133 sec / 2.16888 min / 19.7641 musec/read
	num_mapped_reads 5796427/6584304 (88.034%)

using 8 parallel threads and writing the mapping output to `/dev/null`.

To partition the index to obtain a meta-colored Fulgor index, then do:

	./fulgor meta -i ~/Salmonella_enterica/salmonella_4546.fur -d tmp_dir --check
We can change the first argument from `meta` to `differential` or `meta-differential` to create a differential-colored or a meta-differential-colored index, respectively.

See the table below for some additional data on the different indexes

| command             | output file             | size (GB) | compression factor |
|:--------------------|:------------------------|-----------|--------------------|
| `meta`              | `salmonella_4546.mfur`  | 0.11769   | 2.26               |
| `differential`      | `salmonella_4546.dfur`  | 0.11076   | 2.40               |
| `meta-differential` | `salmonella_4546.mdfur` | 0.09389   | 2.84               |

Index sizes for tested datasets 
-------------------------------
| Dataset |   |       |   | _Fulgor_     |       |   |  _d-Fulgor_  |       |   |  _m-Fulgor_ |       |   | _md-Fulgor_ |       |
|---------|---|-------|---|--------------|-------|---|:------------:|-------|---|:-----------:|-------|---|-------------|-------|
|         |   | dBG   |   | Color sets   | Total |   | Color sets   | Total |   | Color sets  | Total |   | Color sets  | Total |
| EC      |   | 0.29  |   | 1.36   (83%) | 1.65  |   | 0.45   (61%) | 0.74  |   | 0.40  (58%) | 0.69  |   | 0.24  (45%) | 0.52  |
| SE-5k   |   | 0.16  |   | 0.59   (79%) | 0.75  |   | 0.20   (56%) | 0.36  |   | 0.16  (50%) | 0.32  |   | 0.11  (40%) | 0.27  |
| SE-10k  |   | 0.35  |   | 1.66   (83%) | 2.01  |   | 0.48   (58%) | 0.83  |   | 0.34  (49%) | 0.70  |   | 0.22  (39%) | 0.57  |
| SE-50k  |   | 1.25  |   | 17.03  (93%) | 18.29 |   | 4.31   (77%) | 5.57  |   | 2.08  (62%) | 3.34  |   | 1.38  (52%) | 2.64  |
| SE-100k |   | 1.71  |   | 40.71  (96%) | 42.43 |   | 9.37   (84%) | 11.10 |   | 3.75  (68%) | 5.47  |   | 2.26  (57%) | 3.98  |
| SE-150k |   | 2.02  |   | 68.61  (97%) | 70.65 |   | 15.73  (89%) | 17.77 |   | 5.27  (72%) | 7.31  |   | 3.22  (61%) | 5.26  |
| GB      |   | 21.29 |   | 15.54  (42%) | 36.83 |   | 7.51   (26%) | 28.81 |   | 9.16  (30%) | 30.46 |   | 6.19  (23%) | 27.48 |

Pseudoalignment output format
-----------------------------

The tool `pseudoalign` writes the result to an output file, in plain text format, specified with the option `-o [output-filename]`.

This file has one line for each mapped read, formatted as follows:

	[read-name][TAB][list-lenght][TAB][list]

where `[list]` is a TAB-separated list of increasing integers, of length `[list-length]`, representing the list of reference identifiers to which the read is mapped. (`[TAB]` is the character `\t`.)

#### Example

	NODE_11_length_149361_cov_9.71634_ID_21 1       0
	NODE_3406_length_341_cov_20.0437_ID_681	1       0
	NODE_4745_length_118_cov_12.7931_ID_949	3       0       3       7
	NODE_102_length_2047_cov_18.1471_ID_203 1       0
	NODE_477_length_1163_cov_22.0531_ID_953 2       0       8
	NODE_9_length_173161_cov_9.33695_ID_17  1       0
	NODE_22_length_45757_cov_12.1361_ID_43  1       0

#### Important note

If pseudoalignment is performed against a **meta-colored** Fulgor index,
the reference identifiers in the pseudoalignment output might **not** correspond to the ones assigned following the input-file order as specified with option `-l` during index construction.
This is because the meta-colored index re-assignes identifiers to references to improve index compression.

In this case, the reference identifiers in the pseudoalignment output
are consistent with the ones returned by the `print-filenames` tool.