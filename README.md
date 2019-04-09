# transXpress-nextflow
transXpress: a [Nextflow](https://www.nextflow.io) pipeline for rapid de novo transcriptome assembly and annotation

Also see our sister project: [transXpress-snakemake](https://github.com/transXpress/transXpress-snakemake)

## Intro

## Dependencies

Requires
* NextFlow 19.01.0+ (install via conda)
* fastqc (install via conda)
* trimmomatic (install via conda)
* Trinity (install via conda)
* SPAdes (install via conda)
* TransDecoder (install via conda)
* BioPython (install via conda)
* samtools (install via conda)
* bowtie2 (install via conda)
* infernal (install via conda)
* HMMER (install via conda)
* kallisto (install via conda)
* NCBI BLAST+ (install via conda)
* R (install via conda)
* seqkit (install via conda)
* [deeploc](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?deeploc) / [SignalP 4.1](http://www.cbs.dtu.dk/cgi-bin/sw_request?signalp+4.1) / [SignalP 5.0](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?signalp)
* [tmhmm v. 2.0](http://www.cbs.dtu.dk/services/TMHMM/)
* basic Linux utitilies: wget, split, awk, cut, gzip

## Installation

1. Install [Miniconda3](https://conda.io/en/latest/miniconda.html)
2. Install conda dependencies:  
~~~
 conda config --add channels bioconda
 conda config --add channels conda-forge
 conda config --add channels r
 conda install nextflow fastqc trimmomatic trinity=2.8.4-0 spades transdecoder biopython samtools bowtie2 infernal hmmer kallisto blast r seqkit bioconductor-edger
~~~
(Note, below dependencies are optional, transXpress will run to completion without them, but will produce empty files for their output)
3. Install deeploc (performance being evaluated by transXpress developers in comparison to SingalP 4.1/5.0)
      * Download deeploc from http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?deeploc
      * Install dependencies: `pip install -r requirements.txt`
      * Install deeploc: `python setup.py install` or locally: `python setup.py install --user`
      
4. Install SignalP 4.1g (performance being evaluated by transXpress developers in comparison to SingalP 5.0/deeploc)
      * Download SignalP 4.1g from http://www.cbs.dtu.dk/cgi-bin/sw_request?signalp+4.1


5. Install SignalP 5.0 (performance being evaluated by transXpress developers in comparison to SingalP 4.1/deeploc)
      * Download SignalP 5.0 from http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?signalp
      
6. Install tmhmm
      * Download tmhmm from http://www.cbs.dtu.dk/services/TMHMM/


## Usage
Make your assembly directory and change it to the current directory
```
mkdir your_assembly_directory
cd your_assembly_directory
```
Setup the mandatory 'samples.txt' file in the assembly directory describing where to find your raw read FASTQ files. Reads will be pooled from all samples for a single transcriptome assembly, but expression quantification will happen on a per-sample basis. See the tests directory for an example of a samples file: [samples.txt](./tests/test_nonSS-trinity/samples.txt)

Setup the mandatory 'species.txt' file in the directory describing which species the data comes from. See the tests directory for an example of a species file: [species.txt](./tests/test_nonSS-trinity/species.txt)

Symbolically link the transxpress-nextflow code into your assembly directory
```
ln -s /your/transxpress-nextflow-cloned-directory/* ./
```
Execute the run.sh script with your assembler of choice, either `trinity` or `rnaspades` currently
```
./run.sh trinity
```
NextFlow only likes 1 assembly per directory, so if you'd like to run two assemblies simultaneously, you have to use different assembly directories.

## Flow graph
![Directed acyclic graph for transXpress-nextflow program execution](./tests/test_nonSS-trinity/test_nonSS_dag.svg)

