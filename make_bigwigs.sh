#!/bin/bash

#SBATCH --job-name=bigwig
#SBATCH -A ecoevo283                  ## account to charge
#SBATCH -p standard                   ## standard partition
#SBATCH -N 1                          ## nodes
#SBATCH --array=1-24                    ## tasks
#SBATCH -c 2                          ## CPUs
#SBATCH -t 1-                         ## day run time limit


##############################
## DO BEFORE RUNNING SCRIPT ##
##############################
#mkdir analysis
#ls ../alignments/*RG* | cut -d "/" -f 3 > atac_bams_list.txt
#cat $ref.fai | head -n 7 | awk '{printf("%s\t0\t%s\n",$1,$2)}' > major.bed

module load samtools/1.10
module load ucsc-tools/v393 
module load bedtools2/2.29.2

ref="/data/class/ecoevo283/tzhuravl/RAWDATA/ref/dmel-all-chromosome-r6.13.fasta"
file="/data/class/ecoevo283/tzhuravl/RAWDATA/ATACseq/analysis/atac_bams_list.txt"
bam=`head -n $SLURM_ARRAY_TASK_ID  ${file} | tail -n 1`

path="/data/class/ecoevo283/tzhuravl/RAWDATA/ATACseq/analysis"

cd $path


samtools view -b -L major.bed ../alignments/$bam > $bam.maj

Nreads=`samtools view -c -F 4 $bam.maj`

Scale=`echo "1.0/($Nreads/1000000)" | bc -l`

samtools view -b $bam.maj | genomeCoverageBed -bg -scale $Scale -ibam - > $bam.coverage

bedSort $bam.coverage $bam.sort.coverage

bedGraphToBigWig $bam.sort.coverage $ref.fai $bam.bw
