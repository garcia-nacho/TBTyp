# TBTyp   
*Mycobacterium tuberculosis* classification pipeline.   
[![Linux](https://svgshare.com/i/Zhy.svg)](https://svgshare.com/i/Zhy.svg)   [![Docker](https://badgen.net/badge/icon/docker?icon=docker&label)](https://https://docker.com/)
   
This pipeline can be used to classify M. tuberculosis samples from Illumina paired-end fastq files using the [TBProfiler](https://github.com/jodyphelan/TBProfiler) scheme.   
If there is a mixture of lineages based on the different SNPs of the scheme the ratio of reads supporting each of the lineages is displayed.    
   
## Running TBTyp
To run TBTyp you just need docker installed and run this command inside a master folder that contains one folder with two fastq.gz files per read.   

<code>docker run -it --rm -v $(pwd):/Data ghcr.io/garcia-nacho/tbtyp</code>   
   
Expected folder structure:

<pre>
./Experiment         
  |-Sample1     
      |-Sample1_R1.fastq.gz       
      |-Sample1_R2.fastq.gz

  |-Sample2      
      |-Sample2_R1.fastq.gz       
      |-Sample2_R2.fastq.gz
  |-...
  |-...
      
</pre>

## Results   
TBTyp will produce an excel file with the predicted lineage and lineage probabilities. The average depth and certainty (1-Noise) of the SNPs ananlyzed are shown   
It will produce one pdf for sample showing the SNPs found, the lineage assigned to each of them, the depth and the certainty to support the call.


