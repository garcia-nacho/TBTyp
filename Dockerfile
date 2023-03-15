FROM garcianacho/fhibase:v1
LABEL maintainer="Nacho Garcia <iggl@fhi.no>"

USER docker
RUN cd /home/docker \
    && wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /home/docker/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh \
    && conda config --add channels defaults \
    && conda config --add channels bioconda \
    && conda config --add channels conda-forge \
    && conda install -c bioconda samtools\
    && conda install -c bioconda seqkit \
    && conda update -n base -c defaults conda \
    && conda install -c bioconda bedtools \
    && conda install -c bioconda bowtie2 \
    && ln -s /home/docker/miniconda3/lib/libcrypto.so.1.1 /home/docker/miniconda3/lib/libcrypto.so.1.0.0 

USER root 
RUN Rscript -e "install.packages(c('ggplot2','writexl', 'ggrepel', 'data.table','seqinr'))"
RUN mkdir -p /Data /home/docker/CommonFiles
COPY CommonFiles/ /home/docker/CommonFiles/
COPY bin/FINex /usr/bin/FINex
RUN chmod -R +rwx /home/docker/CommonFiles/* \
    && chmod 777 /Data && chmod +x /usr/bin/FINex
USER docker
WORKDIR /Data
CMD ["sh", "-c", "/home/docker/CommonFiles/Scripts/tbtyp.sh"]


