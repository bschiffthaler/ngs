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
RUN apt update
RUN apt install -y \
    wget libghc-zlib-dev \
    libncurses-dev git unzip ftp libzmq3-dev \
    nano python2.7-dev ftp \
    apache2 openssl curl libboost-all-dev \
    supervisor openssh-server libbz2-dev \
    liblzma-dev automake cmake && \
    mkdir build && \
    usermod -G training,www-data training
WORKDIR /build

#########
### HTSlib tools: Samtools, Bcftools
#########
RUN wget https://github.com/samtools/samtools/releases/download/1.4/samtools-1.4.tar.bz2 \
    https://github.com/samtools/bcftools/releases/download/1.4/bcftools-1.4.tar.bz2 \
    https://github.com/samtools/htslib/releases/download/1.4/htslib-1.4.tar.bz2 && \
    tar -xf htslib-1.4.tar.bz2
WORKDIR htslib-1.4
RUN ./configure && make -j 4 && make install
WORKDIR /build
RUN tar -xf bcftools-1.4.tar.bz2
WORKDIR bcftools-1.4
RUN make -j 4 && make install
WORKDIR /build
RUN tar -xf samtools-1.4.tar.bz2
WORKDIR samtools-1.4
RUN make -j 4 && make install
WORKDIR /build

#########
### STAR
#########
RUN wget https://github.com/alexdobin/STAR/archive/2.5.3a.tar.gz && \
    tar -xf 2.5.3a.tar.gz
WORKDIR STAR-2.5.3a/bin/Linux_x86_64_static
RUN cp STAR /usr/local/bin
WORKDIR /build

#########
### kallisto
#########
RUN wget https://github.com/pachterlab/kallisto/releases/download/v0.43.1/kallisto_linux-v0.43.1.tar.gz && \
    tar -xf kallisto_linux-v0.43.1.tar.gz
WORKDIR kallisto_linux-v0.43.1
RUN cp kallisto /usr/local/bin
WORKDIR /build

#########
### Salmon
#########
WORKDIR /build
RUN wget https://github.com/COMBINE-lab/salmon/archive/v0.8.2.tar.gz && \
    tar -xf v0.8.2.tar.gz && \
    cd /build/salmon-0.8.2 && mkdir build && cd build && \
    cmake .. && make -j 4 && make install && cp ../lib/lib* /usr/local/lib/. && \
    cp ../bin/salmon /usr/local/bin/.

#########
### BEDtools
#########
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.26.0/bedtools-2.26.0.tar.gz && \
    tar -xf bedtools-2.26.0.tar.gz && cd bedtools2 && \
    make -j 4 && make install

#########
### Java tools: Trimmomatic, FastQC
#########
RUN wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip && \
    unzip fastqc_v0.11.5.zip && \
    chmod +x FastQC/fastqc && \
    cp -r FastQC /usr/share/ && \
    ln -s /usr/share/FastQC/fastqc /usr/bin/ && \
    wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip && \
    unzip Trimmomatic-0.36.zip && \
    cp -r Trimmomatic-0.36 /usr/share/ && \
    echo '#!/bin/bash\njava -jar /usr/share/Trimmomatic-0.36/trimmomatic-0.36.jar $@' > /usr/bin/trimmomatic && \
    chmod +x /usr/bin/trimmomatic

#########
### SortmeRNA
#########
RUN wget https://github.com/biocore/sortmerna/archive/2.1.tar.gz && \
    tar -xf 2.1.tar.gz
WORKDIR sortmerna-2.1
RUN ./configure && make -j 4 && make install || true
RUN cp indexdb_rna /usr/local/bin && \
    cp -r rRNA_databases/ /usr/share/
WORKDIR /usr/share/rRNA_databases
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
RUN wget https://github.com/lh3/bwa/releases/download/v0.7.15/bwa-0.7.15.tar.bz2 && \
    tar -xf bwa-0.7.15.tar.bz2
WORKDIR bwa-0.7.15
RUN make -j 4 && cp bwa /usr/local/bin

#########
### Pip installs: HTSeq, PySam, MACS2
#########
RUN pip install numpy && \
    pip install htseq && \
    pip install PySam && \
    pip install MACS2 && \
    pip install tornado==4.4.3
RUN apt install -y libfreetype6-dev libpng-dev libxft-dev && \
    pip install multiqc
WORKDIR /build

#########
### USeq and Sissrs
#########
RUN mkdir /data && \
    wget https://github.com/HuntsmanCancerInstitute/USeq/releases/download/USeq_9.0.8/USeq_9.0.8.zip && \
    unzip USeq_9.0.8.zip && mv USeq_9.0.8 /usr/share/USeq && \
    wget http://dir.nhlbi.nih.gov/papers/lmi/epigenomes/sissrs/sissrs_v1.4.tar.gz && \
    tar -xf sissrs_v1.4.tar.gz && \
    cp sissrs.pl /usr/local/bin

#########
### JBrowse to easily view biological data formats
#########
RUN rm /var/www/html/index.html && \
    ln -sf /home/training /var/www/html/home
WORKDIR /var/www/html
RUN wget http://jbrowse.org/releases/JBrowse-1.12.1.zip && \
    unzip JBrowse-1.12.1.zip
WORKDIR /var/www/html/JBrowse-1.12.1
RUN ./setup.sh && \
    mkdir -p data/bam && \
    chown -R www-data:www-data /var/www/html && \
    chmod a+rwx /var/www/html/JBrowse-1.12.1/data/bam && \
    touch data/tracks.conf && \
    chmod a+rw    /var/www/html/JBrowse-1.12.1/data/tracks.conf
WORKDIR /
RUN rm /var/www/html/JBrowse-1.12.1.zip && \
    rm -rf /build/*
ADD scripts/add_JBrowse_tracks.sh /usr/local/bin/add_JBrowse_tracks.sh

#########
### GateOne SSH interface
#########
RUN git clone https://github.com/liftoff/GateOne/ && \
    cd GateOne && python setup.py install && \
    python run_gateone.py --configure && cd ..

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
EXPOSE 80 8888 8787 443
CMD ["/usr/bin/supervisord","-c","/etc/supervisor/conf.d/supervisord.conf"]
