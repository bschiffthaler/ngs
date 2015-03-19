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
	python-matplotlib
RUN mkdir build
WORKDIR /build
RUN wget https://github.com/samtools/samtools/releases/download/1.2/samtools-1.2.tar.bz2 \
	https://github.com/samtools/bcftools/releases/download/1.2/bcftools-1.2.tar.bz2 \
	https://github.com/samtools/htslib/releases/download/1.2.1/htslib-1.2.1.tar.bz2
RUN tar -xf htslib-1.2.1.tar.bz2
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
RUN wget https://github.com/alexdobin/STAR/archive/STAR_2.4.0j.tar.gz
RUN tar -xf STAR_2.4.0j.tar.gz
WORKDIR STAR-STAR_2.4.0j/bin/Linux_x86_64_static
RUN cp STAR /usr/local/bin
WORKDIR /build
RUN wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.2.zip
RUN unzip fastqc_v0.11.2.zip && \
	chmod +x FastQC/fastqc && \
	cp -r FastQC /usr/share/ && \
	ln -s /usr/share/FastQC/fastqc /usr/bin/
RUN wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.33.zip
RUN unzip Trimmomatic-0.33.zip && \
	cp -r Trimmomatic-0.33 /usr/share/
RUN echo '#!/bin/bash\njava -jar /usr/share/Trimmomatic-0.33/trimmomatic-0.33.jar $@' > /usr/bin/trimmomatic && \
	chmod +x /usr/bin/trimmomatic
RUN wget http://bioinfo.lifl.fr/RNA/sortmerna/code/sortmerna-2.0-linux-64.tar.gz && \
	tar -xf sortmerna-2.0-linux-64.tar.gz
WORKDIR sortmerna-2.0-linux-64
RUN cp sortmerna /usr/local/bin && cp indexdb_rna /usr/local/bin
RUN cp -r rRNA_databases/ /usr/share/
WORKDIR /usr/share/rRNA_databases
COPY indexdb.sh /usr/share/rRNA_databases/
RUN ./indexdb.sh
WORKDIR /build
RUN wget http://downloads.sourceforge.net/project/bio-bwa/bwa-0.7.12.tar.bz2 && \
	tar -xf bwa-0.7.12.tar.bz2
WORKDIR bwa-0.7.12
RUN make && \
	cp bwa /usr/local/bin/
RUN pip install htseq
RUN pip install PySam
WORKDIR /
#RUN rm -rf /build
RUN mkdir /data
EXPOSE 8787
CMD ["/usr/bin/supervisord","-c","/etc/supervisor/conf.d/supervisord.conf"]
