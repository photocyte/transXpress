# transXpress-nextflow
transXpress: a Nextflow pipeline for rapid de novo transcriptome assembly and annotation

Also see our sister project: [transXpress-snakemake](https://github.com/transXpress/transXpress-snakemake)

## Intro

## Dependencies

Requires
* NextFlow 19.02.0+
* BioPython
* samtools
* R
* infernal
* seqkit
* basic linux utitilies: wget, split

## Installation


## Usage
Make your assembly directory and change it to the current directory
```
mkdir your_assembly_directory
cd your_assembly_directory
```
Setup the mandatory 'samples.txt' file in the assembly directory describing where to find your raw read FASTQ files. See the tests directory for an example:[samples.txt](./tests/test_nonSS-trinity/samples.txt)

Setup the mandatory 'species.txt' file in the directory describing which species the data comes from. See the tests directory for an example:[species.txt](./tests/test_nonSS-trinity/species.txt)

Symbolically link the transxpress-nextflow code into your assembly directory
```
ln -s /your/transxpress-nextflow-cloned-directory/* ./
```
Execute the run.sh script with your assembler of choice, either trinity or rnaspades currently
```
./run.sh trinity
```
NextFlow only likes 1 assembly per directory, so if you'd like to run two assemblies simultaneously, you have to use different assembly directories.

## Flow

**Trinity**
![Directed acyclic graph for Trinity transXpress-nextflow program execution](./tests/test_nonSS-trinity/test_nonSS_dag.svg)

**rnaSPAdes**
![Directed acyclic graph for rnaSPAdes transXpress-rnaspades program execution](./tests/test_nonSS-rnaspades/test_nonSS_dag.svg)
