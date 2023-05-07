#!/usr/bin/env bash

################################################
# PROGRAM:
# RNAseq_analyzer_221126.sh
#
# DESCRIPTION:
# This is a Tuxedo based RNA-seq pipeline that I wrote for analyzing fq reads. Step1 is a
# simple wrapper that performs quality control, genome alignment, basic format
# conversions, and htseq-count tabulation for paired-end RNA-seq samples using a specified
# genome. This pipeline has also been written to generate a standard ballgown input directory which can be imported into R/ballgown for differential #expression or other downstream analysis. Step2, which is still being written, is a clean up program that removes unnecessary files and compresses #files to save space.
#
# AUTHOR:
# Daniel Kunk (Colorado State University)
#
# START DATE:
# April 4th, 2023
#
# DEPENDENCIES:
# 	Requires the installation of the follwing software: 
#		fastp
#		hisat2
#		Stringtie
#		samtools
#		bedtools
#		deep-tools
#
# 	Requires access to the Nishimura Lab installed software sources on SUMMIT.
#
# REQUIRES:
#    INPUT: .fastq files.    For each sample, paired forward and reverse sequencing files
#								are required. These should be placed in an input
#								directory.
#
#    INPUT: _metadata.txt file: A metadata file with two columns. The first two columns
#								are fastq file names. The third column is a "nickname"
#								of each sample. Later columns can be included with other
#								metadata information. Metadata file should be placed
#								within the inputdir directory.
#
#
#    HISAT2 INDEXES: .ht2 files for the genome. These are produced using hisat2-build. For
#								instructions see
#	           https://ccb.jhu.edu/software/hisat2/manual.shtml#the-hisat2-build-indexer
#
#    GENOME SEQUENCE: .fa  or .tar.gz file for the genome. This is the sequence of the 
#                                genome.
#
#    GENOME ANNOTATION: .gtf file for the genome. This is a genome annotation file of gene
#								features. Version and coordinates must match the genome
#								sequence (.fa above).
#
# USAGE:
# $ bash RNAseq_analyzer_221126.sh <metadata.txt> <number of threads>
#
# OUTPUT:
#
# KNOWN BUGS:
#
# THINGS TO IMPROVE:
#
################################################


####### MODIFY THIS SECTION #############

#The input samples live in directory:
inputdir="/scratch/alpine/dkunk1@colostate.edu/RNAseq_ana_2023/samples/"

#Metadata file. This pulls the metadata path and file from the command line
metadata=$1

#This is where the ht2 files live:
hisat2path="/scratch/alpine/dkunk1@colostate.edu/RNAseq_ana_2023/indexes/rpv2"

#This is where the genome sequence lives:
genomefa="/scratch/alpine/dkunk1@colostate.edu/RNAseq_ana_2023/genome/R_padi_v2.fasta"

#This is where the gtf file lives:
gtffile="/scratch/alpine/dkunk1@colostate.edu/RNAseq_ana_2023/genome/Rpadi_v2.gff"

#This is the output_directory:
DATE=`date +%Y-%m-%d`
#OR
#DATE='2022-12-03'

outputdir="/scratch/alpine/dkunk1@colostate.edu/RNAseq_ana_2023/output_newTuxedo/"


########## DONE MODIFYING ###############



########## BEGIN CODE ###############

echo -e ">>> INITIATING analyzer with command:\n\t$0 $@"

#Number of threads to use:
# Note - this imports the number of threads (ntasks) given in the command line
pthread=$2

# Make output directories
echo -e ">>> MAKING output directory"
echo -e "\tmkdir $outputdir"
mkdir -p $outputdir



####### META DATA #############

#These are the sample names, R1:
samples1=( $(cut -f 1 --output-delimiter=' ' $metadata) )

#These are the sample names, R2:
samples2=( $(cut -f 2 --output-delimiter=' ' $metadata) )

#These are the nicknames I want to give the files:
names=( $(cut -f 3 --output-delimiter=' ' $metadata) )



####### PIPELINE ##############

# Report back to the user which files will be processed and which names they'll be given:
echo -e ">>> INPUT: This script will process files from the metafile:\n\t$metadata"
echo -e ">>> PLAN: This script will process the sample files into the following names: "
echo -e "\tSAMPLE1\tSAMPLE2\tNAMES"

for (( counter=0; counter < ${#samples1[@]}; counter++ ))
do
    echo -e "\t${samples1[$counter]}\t${samples2[$counter]}\t${names[$counter]}"
done


# FASTP to remove unwanted sequences
# FASTP to determine quality
echo -e "\n>>> FASTP: Trimming excess and low-quality sequences from .fastq file; generating quality report"
mkdir -p $outputdir"01_fastp"

for (( counter=0; counter < ${#samples1[@]}; counter++ ))
do
    samplename=${names[$counter]}
    sample1=${samples1[$counter]}
    sample2=${samples2[$counter]}

    ## Echo statements
    
            ##### ENTER ECHO STATEMENTS HERE #####
    
    ## Make output directories
    mkdir -p $outputdir"01_fastp/"$samplename
    
    ## execute fastp
    cmd1="fastp -i $inputdir/$sample1 \
-I $inputdir/$sample2 \
-o ${outputdir}01_fastp/${samplename}/${samplename}_trim_1.fastq \
-O ${outputdir}01_fastp/${samplename}/${samplename}_trim_2.fastq \
-h ${outputdir}01_fastp/${samplename}/${samplename}_report.html \
-j ${outputdir}01_fastp/${samplename}/${samplename}_report.json \
--detect_adapter_for_pe \
--thread $pthread \
-x -g "
    
    echo -e "\t$ ${cmd1}"
    time eval $cmd1

done

# HISAT2 to align to the genome
echo -e "\n>>> HISAT2: aligning each sample to the genome"
outhisat2=$outputdir"02_hisat2/"
mkdir -p $outhisat2

for (( counter=0; counter < ${#samples1[@]}; counter++ ))
do
    samplename=${names[$counter]}
    sample1=${samples1[$counter]}
    sample2=${samples2[$counter]}


    ## execute hisat2
    cmd3="hisat2 -x $hisat2path -1 $outputdir"01_fastp/"$samplename/$samplename"_trim_1.fastq" -2 $outputdir"01_fastp/"$samplename/$samplename"_trim_2.fastq" -S ${outhisat2}${samplename}.sam --summary-file ${outhisat2}${samplename}_summary.txt --no-unal -p $pthread"
    echo -e "\t$ $cmd3"
    time eval $cmd3

done



# SAMTOOLS and BAMCOVERAGE: to convert .sam output to uploadable .bam and .wg files
echo -e "\n>>> SAMTOOLS/BAMCOVERAGE: to convert files to uploadable _sort.bam and _sort.bam.bai files:"
samout=$outputdir"03_samtools/"
mkdir -p $samout

for seqname in ${names[@]}
do
    # echo
    echo -e "\tSamtools and BamCoverage convert: ${seqname}"
    
    # Samtools: compress .sam -> .bam
    cmd4="samtools view --threads $pthread -bS ${outhisat2}${seqname}.sam > ${samout}${seqname}.bam"
    echo -e "\t$ ${cmd4}"
    time eval $cmd4

    
    # Samtools: sort .bam -> _sort.bam
    cmd5="samtools sort --threads $pthread -o ${samout}${seqname}_sort.bam --reference $genomefa ${samout}${seqname}.bam"
    echo -e "\t$ ${cmd5}"
    time eval $cmd5
    
    
    # Samtools: index _sort.bam -> _sort.bam.bai
    cmd6="samtools index ${samout}${seqname}_sort.bam"
    echo -e "\t$ ${cmd6}"
    time eval $cmd6
    
    
    # bamCoverage: Create a .bw file that is normalized. This can be uploaded to IGV or UCSC
    cmd7="bamCoverage -b ${samout}${seqname}_sort.bam -o ${samout}${seqname}_sort.bw --outFileFormat bigwig -p $pthread --normalizeUsing CPM --binSize 1"
    echo -e "\t$ ${cmd7}"
    time eval $cmd7
    
done



# String-tie: To assemble and quantify expressed genes and transcripts:
 
echo -e "\n>>> String-tie: To assemble and quantify expressed genes and transcripts:"
stringtieout=$outputdir"04_string-tie/"
mkdir -p $stringtieout

for seqname in ${names[@]}
do  
    cmd8="stringtie -p $pthread -G $gtffile -o ${stringtieout}${seqname}.gtf $outputdir"03_samtools/"$seqname"_sort.bam""
    
    echo -e "\t$ ${cmd8}"
    time eval $cmd8
done


# String-tie: To merge all gtf/gff files given by String-tie into one uniformn, non-redundant transcriptome. 
# Here we do not have to use a loop since stringtie will take a list of files and generate one merged file.
echo -e "making the gff_files.txt file which is a list of all the gff's"
touch gff_files.txt

# finding all outputted gff files and making a text file with their path and file name and adding them to the new gtf_list.txt file
echo -e "grabbing all the gff file names and putting them in the text document"
find $stringtieout -type f -name "*.gtf" -printf "%p\n" > gff_files.txt

cmd9="stringtie --merge -G $gtffile -p $pthread -o ${stringtieout}Rpadi_head_stringtie_merged.gtf gff_files.txt"

echo -e "\t$ ${cmd9}"
    time eval $cmd9

# String-tie: Now we are trying to create "Ballgown Input Table Files".
# These are tables which contains coverage data for all transcripts.
# It seems that running this with -b and -e will still make the output of -o which is a gtf file (we already made and merged these) so we will just save these new gtf files into a new directory which will be called 05_ballgown

echo -e "\n>>> String-tie: To create Ballgown Input Table Files:"
ballgown=$outputdir"05_ballgown/"
mkdir -p $ballgown

merged=${stringtieout}Rpadi_head_stringtie_merged.gtf

for seqname in ${names[@]}
do
    
    cmd10="stringtie -e $outputdir"03_samtools/"$seqname"_sort.bam" -B -p $pthread -G $merged -o $outputdir"05_ballgown/"$seqname/$seqname".gtf" -A $outputdir"05_ballgown/"$seqname/$seqname"gene_abund.tab"  "
    
    echo -e "\t$ ${cmd10}"
    time eval $cmd10
done

echo -e "RNASEQ ANALYSIS COMPLETE - New Tuxedo Pipeline DK2023"

######## VERSIONS #############
echo -e "\n>>> VERSIONS:"
echo -e "\n>>> FASTP VERSION:"
fastp --version
echo -e "\n>>> HISAT2 VERSION:"
hisat2 --version
echo -e "\n>>> SAMTOOLS VERSION:"
samtools --version
echo -e "\n>>> STRINGTIE VERSION:"
stringtie --version
echo -e "\n>>> BAMCOVERAGE VERSION:"
bamCoverage --version
echo -e ">>> END: Analayzer complete."