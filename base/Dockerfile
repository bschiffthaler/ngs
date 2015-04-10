FROM bioconductor/release_sequencing
MAINTAINER Bastian Schiffthaler (bastian@bioinformatics.upsc.se)
RUN useradd -m training
RUN echo 'training:training' | chpasswd
RUN mkdir /home/training/reference
RUN mkdir /home/training/data
RUN chsh -s /bin/bash training
RUN echo "alias ll='ls -la -G'" >> /home/training/.profile
RUN apt-get update &&  apt-get install -y \
	build-essential wget libghc-zlib-dev \
	libncurses-dev git unzip openjdk-7-jdk \
	nano python2.7-dev python-numpy \
	python-matplotlib apache2 && \
	mkdir build && \
	usermod -G training,www-data training
WORKDIR /build
RUN wget https://github.com/samtools/samtools/releases/download/1.2/samtools-1.2.tar.bz2 \
	https://github.com/samtools/bcftools/releases/download/1.2/bcftools-1.2.tar.bz2 \
	https://github.com/samtools/htslib/releases/download/1.2.1/htslib-1.2.1.tar.bz2 && \
	tar -xf htslib-1.2.1.tar.bz2
WORKDIR htslib-1.2.1
RUN ./configure && make && make install
WORKDIR /build
RUN tar -xf bcftools-1.2.tar.bz2
WORKDIR bcftools-1.2
RUN make && make install
WORKDIR /build
RUN tar -xf samtools-1.2.tar.bz2
WORKDIR samtools-1.2
RUN make && make install
WORKDIR /build
RUN wget https://github.com/alexdobin/STAR/archive/STAR_2.4.0j.tar.gz && \
	tar -xf STAR_2.4.0j.tar.gz
WORKDIR STAR-STAR_2.4.0j/bin/Linux_x86_64_static
RUN cp STAR /usr/local/bin
WORKDIR /build
RUN wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.2.zip && \
	unzip fastqc_v0.11.2.zip && \
	chmod +x FastQC/fastqc && \
	cp -r FastQC /usr/share/ && \
	ln -s /usr/share/FastQC/fastqc /usr/bin/ && \
	wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.33.zip && \
	unzip Trimmomatic-0.33.zip && \
	cp -r Trimmomatic-0.33 /usr/share/ && \
	echo '#!/bin/bash\njava -jar /usr/share/Trimmomatic-0.33/trimmomatic-0.33.jar $@' > /usr/bin/trimmomatic && \
	chmod +x /usr/bin/trimmomatic && \
	wget http://bioinfo.lifl.fr/RNA/sortmerna/code/sortmerna-2.0-linux-64.tar.gz && \
	tar -xf sortmerna-2.0-linux-64.tar.gz
WORKDIR sortmerna-2.0-linux-64
RUN cp sortmerna /usr/local/bin && cp indexdb_rna /usr/local/bin && \
	cp -r rRNA_databases/ /usr/share/
WORKDIR /usr/share/rRNA_databases
COPY scripts/indexdb.sh /usr/share/rRNA_databases/
COPY scripts/makeDbList.sh /usr/share/rRNA_databases/
RUN ./indexdb.sh
ENV SORTMERNA_DB="/usr/share/rRNA_databases/rfam-5s-database-id98.fasta,/usr/share/rRNA_databases/index/rfam-5s-database-id98:/usr/share/rRNA_databases/silva-euk-28s-id98.fasta,/usr/share/rRNA_databases/index/silva-euk-28s-id98:/usr/share/rRNA_databases/silva-euk-18s-id95.fasta,/usr/share/rRNA_databases/index/silva-euk-18s-id95:/usr/share/rRNA_databases/silva-arc-16s-id95.fasta,/usr/share/rRNA_databases/index/silva-arc-16s-id95:/usr/share/rRNA_databases/silva-bac-16s-id90.fasta,/usr/share/rRNA_databases/index/silva-bac-16s-id90:/usr/share/rRNA_databases/silva-arc-23s-id98.fasta,/usr/share/rRNA_databases/index/silva-arc-23s-id98:/usr/share/rRNA_databases/silva-bac-23s-id98.fasta,/usr/share/rRNA_databases/index/silva-bac-23s-id98:/usr/share/rRNA_databases/rfam-5.8s-database-id98.fasta,/usr/share/rRNA_databases/index/rfam-5.8s-database-id98"
ADD scripts/merge-paired-reads.sh /usr/local/bin/
ADD scripts/unmerge-paired-reads.sh /usr/local/bin/
WORKDIR /build
RUN wget http://downloads.sourceforge.net/project/bio-bwa/bwa-0.7.12.tar.bz2 && \
	tar -xf bwa-0.7.12.tar.bz2
WORKDIR bwa-0.7.12
RUN make && \
	cp bwa /usr/local/bin/ && \
	pip install htseq && \
	pip install PySam
WORKDIR /
RUN mkdir /data
RUN rm /var/www/html/index.html && \
    ln -sf /home/training /var/www/html/home
WORKDIR /var/www/html
RUN curl -O http://jbrowse.org/releases/JBrowse-1.11.6.zip && \
    unzip JBrowse-1.11.6.zip
WORKDIR /var/www/html/JBrowse-1.11.6
RUN ./setup.sh && \
    mkdir -p data/bam && \
    chown -R www-data:www-data /var/www/html && \
    chmod a+rwx /var/www/html/JBrowse-1.11.6/data/bam && \
    touch data/tracks.conf && \
    chmod a+rw	/var/www/html/JBrowse-1.11.6/data/tracks.conf
WORKDIR /
RUN rm /var/www/html/JBrowse-1.11.6.zip && \
    rm -rf /build/*
ADD scripts/add_JBrowse_tracks.sh /usr/local/bin/add_JBrowse_tracks.sh
ADD supervisord.conf /etc/supervisor/conf.d/supervisord.conf
EXPOSE 8787
EXPOSE 80
CMD ["/usr/bin/supervisord","-c","/etc/supervisor/conf.d/supervisord.conf"]
