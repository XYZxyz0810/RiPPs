Bioinformatics Pipeline README
Overview
This is a Python-based bioinformatics pipeline for processing nucleotide sequences, predicting ORFs, performing functional annotation, and clustering short peptides. The tool integrates several bioinformatics programs including Prodigal, HMMER (hmmscan), and MMseqs2.

Features
Fetches nucleotide sequences from NCBI using accession numbers or processes existing files

Predicts open reading frames (ORFs) using Prodigal

Annotates protein domains using HMMER against Pfam database

Clusters short peptides using MMseqs2

Generates sequence logos for peptide clusters

Processes cblaster output files

Requirements
Python 3.x

Linux system

Required Python packages (install with pip install -r requirements.txt):

argparse

os

shutil

External Dependencies
The following external programs must be installed and available in your PATH:

Prodigal - ORF prediction

HMMER - protein domain annotation (specifically hmmscan)

MMseqs2 - peptide clustering

Installation
Clone this repository:

bash
git clone [your-repository-url]
cd [repository-name]
Download the Pfam-A.hmm database file and place it in the main directory or specify its path when running.

Usage
The pipeline can be run in three different modes depending on your input:

1. Using an accession file
bash
python3 Main.py --accession_file accession_file_example.txt --updown_distance 10000 --pre_min 20 --pre_max 150 --out_dir output_directory
2. Using an existing nucleotide file
bash
python3 Main.py --nucl_file nucl_file_example.txt --updown_distance 10000 --pre_min 20 --pre_max 150 --out_dir output_directory
3. Using a cblaster result file
bash
python3 Main.py --cblaster_file cblaster_file_example.txt --updown_distance 10000 --pre_min 20 --pre_max 150 --out_dir output_directory
Arguments
Argument	Description
-a, --accession_file	Input file containing NCBI accession numbers
-ud, --updown_distance	Upstream and downstream DNA length to fetch (required)
--pre_min, -pmin	Minimal length of short peptide (integer)
--pre_max, -pmax	Maximal length of short peptide (integer)
--out_dir, -o	Output directory for results
--nucl_file, -n	Path to an existing nucleotide file
--cblaster_file, -c	Binary result file from cblaster
Output Files
The pipeline generates several output files in the specified output directory:

Nucleotide sequence files (.fasta)

Prodigal prediction files (.faa, .gbk)

CDS and peptide prediction files

HMMER domain annotation results

Cluster prediction tables

MMseqs2 clustering results

Sequence logos for peptide clusters

Example Files
Example input files are provided in the repository:

accession_file_example.txt - Example accession number file

nucl_file_example.txt - Example nucleotide sequence file

cblaster_file_example.txt - Example cblaster result file
