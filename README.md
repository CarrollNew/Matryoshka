Matryoshka

Matryoshka is a software for prediction transposon structure in bacterial genomes. Our pipeline based on prediction of transposase genes via Hidden Markov models (HMMs) and finding the repeats around it.


1. Quick instructions on installing Matryoshka:

1. Install Prodigal (https://github.com/hyattpd/Prodigal)
2. Install HMMER (http://hmmer.org/)
3. Install Matryoshka via this command:


    git clone https://github.com/CarrollNew/Matryoshka


2. Running of Matryoshka

For running Matryoshka you should upload the reference genome of bacteria you interested in fasta format to main directory and run script via this command:


    bash everything.sh <reference_bacterial_genome.fasta>

After this procedure program will generate output files:


- reference_bacterial_genome.fasta.gdf_Information.csv - description of transposons found using Matryoshka. This file contains information about coordinates of 5’ and 3’ flanking repeats at the ends of transposon, length of these repeats and information about genes located inside the transposon. Our program annotates genes as transposases and non-transposases. This file also contains information about gene id’s of non-transposases located in transposon.
- reference_bacterial_genome.fasta.faa - protein sequences of genes located in predicted transposons.
