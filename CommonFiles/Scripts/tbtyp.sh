#!/bin/bash

mkdir /home/docker/results


for dir in $(ls -d */)
do

    cd ${dir}
    echo Analyzing ${dir%/}
    
    R1=$(ls *R1*)
    R2=$(ls *R2*)
    
    (bowtie2 -p 1 -x /home/docker/CommonFiles/Refs/MTBref_bt2 -1 ${R1} -2 ${R2} -S ${dir%/}.sam) 2> ${dir%/}_Bowtie2summary.txt
    samtools view -b -o ${dir%/}.bam  ${dir%/}.sam
    samtools sort ${dir%/}.bam -o ${dir%/}.sorted.bam
    samtools index ${dir%/}.sorted.bam
    samtools depth -a ${dir%/}.sorted.bam > ${dir%/}_depth.tsv
    FINex -f ${dir%/}.sorted.bam > ${dir%/}_noise.tsv
    cp ${dir%/}_noise.tsv /home/docker/results/${dir%/}_noise.tsv
    cp ${dir%/}.sorted.bam /home/docker/results/${dir%/}.sorted.bam
    cp ${dir%/}.sorted.bam.bai /home/docker/results/${dir%/}.sorted.bam.bai
    cp ${dir%/}_Bowtie2summary.txt /home/docker/results/${dir%/}_Bowtie2summary.txt
    rm *.tsv
    rm *.bam
    rm *.sam
    rm *.bai
     
    cd ..

done

cp -R /home/docker/results /Data/TBTypResults
cd /Data/TBTypResults

Rscript /home/docker/CommonFiles/LineageAssignment.R

mkdir bam
mv *.bam bam
mv *.bai bam
mkdir noise
mv *.tsv noise
mkdir plots
mv *.pdf plots
mv *.txt bam
