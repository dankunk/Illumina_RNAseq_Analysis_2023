# Illumina_RNAseq_Analysis_2023
Scripts for running RNA seq analysis on paired end Illumina .fq sequences. These scripts use the Tuxedo suite for de novo transcriptome generation and quantification with Salmon. 


The directory structure for the server or computer where these scripts should run is 

      geneome-
      A folder containing the organisms genome file in fasta format (.fa/.fasta) and the organisms reference file (GTF/GFF).

      Towards the end of the pipeline, this folder will also contain the final transcriptome assembly which will be created with gffread. A index of the initial genome         fasta will be produced during this step. This can be done by running the gff_to_fasta.sbatch script. 

      indexes-
      A folder containing the ht2 indexes built from the organisms genome. 
      These can be create by running the dk_build_aphid_indices.sbatch script in this directory prior to running the pipeline.

      output_newTuxedo- 
      A folder which will be populated during the pipeline running. There will be 01_fastp(trimmed reads), 02_hisat2 (sam and summaries), 03_samtools (bam, bam index,         and bw files), 04_stringtie (each individual gff per sample and one merged gff containing all novel isoforms), 05_ballgown (ballgown output directory for each           sample, even though the rest of the pipeline will be done w/ salmon and the new transcriptome).

      samples-
      A folder containing the PE .fq files which have been pre-unzipped using the unzip_script.sh and the execute_unzip.sbatch scripts. I will try and just use one             script for this operation so that other can more easily complete this task. This script should be run inside of this directory for it to work. It can also be done       in the command line, however the use of more CPU power will speed up the unzip process. 

      scripts- 
      A folder containing the scripts used in the pipeline. You need to be inside of this directory to run most of these scripts. The scripts were written with direct         paths to any files needed so there is no need to bring these into the working directory unless specified otherwise.
