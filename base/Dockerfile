FROM bioconductor/release_core2
MAINTAINER Bastian Schiffthaler (b.schiffthaler@gmail.com)
RUN groupadd -g 2000 training && useradd -m -u 2000 -g 2000 training
RUN echo 'training:training' | chpasswd
RUN mkdir /home/training/reference
RUN mkdir /home/training/data
RUN chsh -s /bin/bash training
RUN echo "alias ll='ls -la -G'" >> /home/training/.profile

#########
### Aptitude packages
#########
RUN apt update && apt install -y --fix-missing \
    wget libghc-zlib-dev libfreetype6-dev libpng-dev libxft-dev \
    libncurses-dev git unzip ftp libzmq3-dev \
    nano python3-dev python3-pip ftp \
    apache2 openssl curl libboost-all-dev \
    supervisor openssh-server libbz2-dev \
    liblzma-dev automake cmake openjdk-8-jre-headless shellinabox && \
    mkdir build && \
    usermod -G training,www-data training
WORKDIR /build

#########
### HTSlib tools: Samtools, Bcftools
#########
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 \
    https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2 \
    https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 && \
    tar -xf htslib-1.9.tar.bz2
WORKDIR htslib-1.9
RUN ./configure && make -j 4 && make install
WORKDIR /build
RUN tar -xf bcftools-1.9.tar.bz2
WORKDIR bcftools-1.9
RUN make -j 4 && make install
WORKDIR /build
RUN tar -xf samtools-1.9.tar.bz2
WORKDIR samtools-1.9
RUN make -j 4 && make install
WORKDIR /build

#########
### STAR
#########
RUN wget https://github.com/alexdobin/STAR/archive/2.7.3a.tar.gz && \
    tar -xf 2.7.3a.tar.gz
WORKDIR STAR-2.7.3a/bin/Linux_x86_64_static
RUN cp STAR* /usr/local/bin
WORKDIR /build

#########
### kallisto
#########
RUN wget https://github.com/pachterlab/kallisto/releases/download/v0.46.1/kallisto_linux-v0.46.1.tar.gz && \
    tar -xf kallisto_linux-v0.46.1.tar.gz
WORKDIR /build/kallisto
RUN cp kallisto /usr/local/bin
WORKDIR /build

#########
### Salmon
#########
WORKDIR /build
RUN wget https://github.com/COMBINE-lab/salmon/releases/download/v1.0.0/salmon-1.0.0_linux_x86_64.tar.gz && \
    tar -xf salmon-1.0.0_linux_x86_64.tar.gz
WORKDIR /build/salmon-latest_linux_x86_64/
RUN cp bin/salmon /usr/local/bin && \
    cp lib/liblzma.so.0 /usr/local/lib && \
    cp lib/libtbb* /usr/local/lib

#########
### BEDtools
#########
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.29.0/bedtools-2.29.0.tar.gz && \
    tar -xf bedtools-2.29.0.tar.gz && cd bedtools2 && \
    make -j 4 && make install

#########
### Java tools: Trimmomatic, FastQC
#########
RUN wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip && \
    unzip fastqc_v0.11.8.zip && \
    chmod +x FastQC/fastqc && \
    cp -r FastQC /usr/share/ && \
    ln -s /usr/share/FastQC/fastqc /usr/bin/ && \
    wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip && \
    unzip Trimmomatic-0.39.zip && \
    cp -r Trimmomatic-0.39 /usr/share/ && \
    echo '#!/bin/bash\njava -jar /usr/share/Trimmomatic-0.39/trimmomatic-0.39.jar $@' > /usr/bin/trimmomatic && \
    chmod +x /usr/bin/trimmomatic

#########
### SortmeRNA
#########
WORKDIR /build
COPY sortmerna-3.0.4-Linux.sh /build/sortmerna-3.0.4-Linux.sh
RUN sh sortmerna-3.0.4-Linux.sh --skip-license --prefix=/usr --exclude-subdir
RUN mkdir /usr/share/rRNA_databases
WORKDIR /usr/share/rRNA_databases
RUN wget https://raw.githubusercontent.com/biocore/sortmerna/master/rRNA_databases/rfam-5.8s-database-id98.fasta && \
    wget https://raw.githubusercontent.com/biocore/sortmerna/master/rRNA_databases/rfam-5s-database-id98.fasta && \
    wget https://raw.githubusercontent.com/biocore/sortmerna/master/rRNA_databases/silva-arc-16s-id95.fasta && \
    wget https://raw.githubusercontent.com/biocore/sortmerna/master/rRNA_databases/silva-arc-23s-id98.fasta && \
    wget https://raw.githubusercontent.com/biocore/sortmerna/master/rRNA_databases/silva-euk-28s-id98.fasta && \
    wget https://raw.githubusercontent.com/biocore/sortmerna/master/rRNA_databases/silva-euk-18s-id95.fasta && \
    wget https://raw.githubusercontent.com/biocore/sortmerna/master/rRNA_databases/silva-bac-16s-id90.fasta && \
    wget https://raw.githubusercontent.com/biocore/sortmerna/master/rRNA_databases/silva-bac-23s-id98.fasta
COPY scripts/indexdb.sh /usr/share/rRNA_databases/
COPY scripts/makeDbList.sh /usr/share/rRNA_databases/
RUN ./indexdb.sh
ENV SORTMERNA_DB="/usr/share/rRNA_databases/rfam-5s-database-id98.fasta,/usr/share/rRNA_databases/index/rfam-5s-database-id98:/usr/share/rRNA_databases/silva-euk-28s-id98.fasta,/usr/share/rRNA_databases/index/silva-euk-28s-id98:/usr/share/rRNA_databases/silva-euk-18s-id95.fasta,/usr/share/rRNA_databases/index/silva-euk-18s-id95:/usr/share/rRNA_databases/silva-arc-16s-id95.fasta,/usr/share/rRNA_databases/index/silva-arc-16s-id95:/usr/share/rRNA_databases/silva-bac-16s-id90.fasta,/usr/share/rRNA_databases/index/silva-bac-16s-id90:/usr/share/rRNA_databases/silva-arc-23s-id98.fasta,/usr/share/rRNA_databases/index/silva-arc-23s-id98:/usr/share/rRNA_databases/silva-bac-23s-id98.fasta,/usr/share/rRNA_databases/index/silva-bac-23s-id98:/usr/share/rRNA_databases/rfam-5.8s-database-id98.fasta,/usr/share/rRNA_databases/index/rfam-5.8s-database-id98"
ADD scripts/merge-paired-reads.sh /usr/local/bin/
ADD scripts/unmerge-paired-reads.sh /usr/local/bin/
WORKDIR /build

#########
### BWA
#########
RUN wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2 && \
    tar -xf bwa-0.7.17.tar.bz2
WORKDIR bwa-0.7.17
RUN make -j 4 && cp bwa /usr/local/bin

#########
### Pip installs: HTSeq, PySam, MACS2
#########
RUN pip3 install numpy && \
    pip3 install htseq && \
    pip3 install PySam && \
    pip3 install MACS2
RUN pip3 install multiqc
WORKDIR /build

#########
### USeq and Sissrs
#########
RUN mkdir /data && \
    wget https://github.com/HuntsmanCancerInstitute/USeq/releases/download/USeq_9.2.3/USeq_9.2.3.zip && \
    unzip USeq_9.2.3.zip && mv USeq_9.2.3 /usr/share/USeq && \
    wget http://dir.nhlbi.nih.gov/papers/lmi/epigenomes/sissrs/sissrs_v1.4.tar.gz && \
    tar -xf sissrs_v1.4.tar.gz && \
    cp sissrs.pl /usr/local/bin

#########
### JBrowse to easily view biological data formats
#########
RUN rm /var/www/html/index.html && \
    ln -sf /home/training /var/www/html/home
WORKDIR /var/www/html
RUN wget https://github.com/GMOD/jbrowse/releases/download/1.16.6-release/JBrowse-1.16.6.zip && \
    unzip JBrowse-1.16.6.zip
WORKDIR /var/www/html/JBrowse-1.16.6
RUN ./setup.sh && \
    mkdir -p data/bam && \
    chown -R www-data:www-data /var/www/html && \
    chmod a+rwx /var/www/html/JBrowse-1.16.6/data/bam && \
    touch data/tracks.conf && \
    chmod a+rw    /var/www/html/JBrowse-1.16.6/data/tracks.conf
WORKDIR /
RUN rm /var/www/html/JBrowse-1.16.6.zip && \
    rm -rf /build/*
ADD scripts/add_JBrowse_tracks.sh /usr/local/bin/add_JBrowse_tracks.sh

#########
### R installs
#########
ENV R_LIBS="/home/training/.r-library"
ADD scripts/R_installs.R /build/R_installs.R
ADD .Renviron /home/training/.Renviron
RUN mkdir /home/training/.r-library && \
    chown -R training:training /home/training/.Renviron && \
    chown -R training:training /home/training/.r-library && \
    chmod 755 /home/training/.r-library && \
    sudo -u training Rscript /build/R_installs.R

ADD supervisord.conf /etc/supervisor/conf.d/supervisord.conf
VOLUME /home/training/share
EXPOSE 80 8888 8787 443
CMD ["/usr/bin/supervisord","-c","/etc/supervisor/conf.d/supervisord.conf"]

