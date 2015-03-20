# General information

I am shamelessly copy/pasting the .MD file from the docker hub page which this source is used to build. 

The dockerhub page is [here](https://registry.hub.docker.com/u/bschiffthaler/ngs/)!

Build the docker with:

```
docker build -t <my_image_name> base
```

----------------------

# General information

Ready-to-work docker for next generation sequence analysis including binaries:

* Sequence data QC [(FastQC)](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* Trimming [(Trimmomatic)](http://www.usadellab.org/cms/?page=trimmomatic) [1]
* rRNA filtering [(SortMeRNA)](http://bioinfo.lifl.fr/RNA/sortmerna/) [2]
* Genome mapping ([STAR](https://github.com/alexdobin/STAR) [3] , [BWA](http://bio-bwa.sourceforge.net/) [4] )
* Feature Summarisation [(HTSeq)](http://www-huber.embl.de/HTSeq/doc/overview.html) [5]
* File manipulation and exploration [(samtools,htslib,bcftools)](http://www.htslib.org/) [10] [11]

For downstream analysis, this docker is based on bioconductor/release_sequencing [6] , which contains all the most commonly used downstream analysis tools implemented in R [7] .

The ":with-data" tagged image contains a set of training data which is commonly used in our RNA-Seq training courses. The data is a sub-sampled set from Robinson, Delhomme et al. 2014 [9]. I also included the _Populus trichocarpa_ genome assembly and annotation [8] as well as a _P. tremula_ in-house draft assembly (not published, contact us for more information).

The source (+Dockerfile) which was used to build this container, is [here](https://github.com/bschiffthaler/ngs) on GitHub! 

# Common use cases

Even if meant for training purposes, this Docker should allow a scientist to carry out an RNA-Seq based differential expression study as described [in this protocol](http://www.epigenesys.eu/images/stories/protocols/pdf/20150303161357_p67.pdf).

For reasons of brevity, I will not be going into great detail (see the protocol), but I will give a short overview of common use cases in differential expression analysis.

## Index

1. Basics
	*  Mounting host directories inside the docker container
1. Technical quality control with FastQC
1. rRNA sorting with SortMeRNA
1. Quality trimming with Trimmomatic
1. Mapping to the reference genome with STAR
	*  Creating a STAR genome
	*  Mapping my reads
1. Feature summarisation using HTSeq
1. Downstream analysis using RStudio
1. References

## Basics

### Mounting host directories inside the docker container

In order to analyse your data from within the docker, you need to mount the directory containing it in the docker container. In this example, I have a FASTQ file called my_sample.fq.gz in my current working directory and would like to analyse it with software from my docker. I therefore mount my current working directory in the `/data` folder within the container. The `--rm` flag tells docker to clean up the process after it's done.:

```
docker run --rm -v $(pwd):/data bschiffthaler/ngs zcat /data/my_sample.fq.gz | head
```

The argument `-v $(pwd):/data` tells docker to mount the current working directory (from the `pwd` command) in the folder called `/data` within the container. Then, the `zcat` command is executed, which decompresses the file and streams the contents to the terminal. We pipe (using `|`) the output to the `head` command so that we only read the first few lines. Of course, this example is almost entirely non practical as you could do all of this without having the docker image.

On Microsoft windows, the `$(pwd)` shortcut would not work to point to the current directory, here you can use `CD`.

## Technical quality control with FastQC

As a small prologue to this section I would like to differentiate between technical and biological QA. As an example, technical QA would show sequencing errors, PCR over-amplification and similar issues, while a sample swap during the RNA extraction would not show up during this part. For the lack of a better term, we call methods to determine the latter type of issues "biological QA". This part is about technical QA.

After you receive your sequencing data, and after every step that modifies it (trimming, rRNA sorting) QA should be done. Here, we will be using Simon Andrews' [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/):

```
#Here I assume that I have a FASTQ file called 202_subset_1.fq.gz in my current working directory
docker run --rm -v $(pwd):/data bschiffthaler/ngs fastqc /data/202_subset_1.fq.gz
```

Note that if you're running docker in Microsoft Windows, the `$(pwd)` shortcut will not work and you will have to replace this with the absolute path of the folder you would like to mount.

This will output the FastQC report as html in the directory. For interpretation of that report, please refer to the manual or the aforementioned protocol.

## rRNA sorting with SortMeRNA

Sorting out rRNA contamination with `sortmerna` can be carried out in a similar manner. The rRNA databases which come with `sortmerna` are in `/usr/share/rRNA_databases` and have been pre-indexed during the build of the docker. The indexes are in `/usr/share/rRNA_databases/index` and are always prefixed with the same name (minus the ".fasta" ending) as the source FASTA file. The correct syntax of the command can be taken from `sortmerna -h` and the [sortmerna manual](http://bioinfo.lifl.fr/RNA/sortmerna/code/SortMeRNA-user-manual-v2.0.pdf).

To list the available databases:

```
docker run --rm bschiffthaler/ngs find /usr/share/rRNA_databases -name "*.fasta"
```

For convenience, I included an environment variable $SORTMERNA_DB, which can be passed to `sortmerna`'s `--ref` argument:

```
docker run --rm bschiffthaler/ngs bash -c 'sortmerna --ref $SORTMERNA_DB <...>'
```
Wrapping the command in `bash -c` here is crucial, since otherwise the environment variable would not be passed correctly.

If I would like to sort my FASTQ file which I got from sequencing:

```
docker run -i -t --rm -v $(pwd):/data bschiffthaler/ngs sortmerna --ref \
/usr/share/rRNA_databases/silva-bac-16s-id90.fasta,/usr/share/rRNA_databases/index/silva-bac-16s-id90 \
--reads /data/202_subset_1.fq \
--other /data/202_subset_1_sorted \
--aligned /data/202_subset_1_rRNA_hits --fastx

```

This outputs two files: 

1. 202_subset_1_sorted.fq - which contains all reads which show **NO** hits to an rRNA
2. 202_subset_1_rRNA_hits.fq - which contains all reads which **DO** show hits to an rRNA

## Quality trimming with Trimmomatic

In case quality trimming (i.e.: dropping low quality reads or adapter sequences) is desired, this docker comes with Trimmomatic. I will again leave the specifics to the [Trimmomatic manual](http://www.usadellab.org/cms/?page=trimmomatic) and refer to the protocol for discussion about trimming. In the example case, I decide to trim my already rRNA-sorted data with Trimmomatic. To do so, I run:

```
docker run -i -t --rm -v $(pwd):/data bschiffthaler/ngs trimmomatic SE -phred64 \
/data/202_subset_1_sorted.fq /data/202_subset_1_sorted_trimmed.fq LEADING:10 \
TRAILING:10 SLIDINGWINDOW:5:20 \
ILLUMINACLIP:/usr/share/Trimmomatic-0.33/adapters/TruSeq2-SE.fa:2:30:10
```

This leaves me with a quality trimmed output file called "202_subset_1_sorted_trimmed.fq", which I will now map to the genome in the next step.

## Mapping to the reference genome with STAR

### Creating a STAR genome

After downloading an appropriate reference genome assembly for my organism, I generate a STAR genome which STAR will use during the mapping step. NOTE: This step might require workstation with a somewhat large amount of RAM (see [STAR manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)). If you do not have that available there are several parameters available which will reduce the memory requirement (at the cost of mapping speed). Those are described in the appendix of the STAR manual. If available, it is highly advised to feed the annotation in the form of a GTF file into STAR's `genomeGenerate`. If the annotation is provided as GFF, the `cufflinks` software suite comes with a tool called `gffread` that is capable of converting the formats:

```
docker run -i -t --rm -v $(pwd):/data bschiffthaler/ngs STAR --runMode \
genomeGenerate --genomeDir /data/Potri_star/ --genomeFastaFiles \
/data/Ptrichocarpa_v3.0_210.fa --sjdbGTFfile /data/Ptrichocarpa_v3.0_210_gene_exons.gtf \
--sjdbOverhang 99
```

The folder "Potri_star" (which I created beforehand) now holds the STAR genome.

### Mapping my reads

As always, refer to the manual for exact details on STAR's command line options. To summarise quickly, I direct STAR to my previously created genome directory and my sequencing read files, I set my maximum desired intron length to 11,000 (as this is the longest intron in the annotation) and finally, I specify that I would like the output file in the BAM format and sorted by coordinate.

```
docker run -i -t --rm -v $(pwd):/data bschiffthaler/ngs STAR --genomeDir \
/data/Potri_star --readFilesIn /data/202_subset_1_sorted_trimmed.fq \
--alignIntronMax 11000 --outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix /data/202_subset_1_sorted_trimmed_STAR
```

## Feature summarisation using HTSeq

The final step before the actual differential expression testing is the feature summaristaion using HTSeq, specifically `htseq-count`. In this case we need to run the executable from within `bash -c` like this: `bash -c 'htseq-count <arguments>'`. The reason is that htseq-count does not write to a file, but writes to stdout, which means that we need to redirect stdout to a file. We could do `docker run <...> htseq-count <...> > my_counts.txt`, but only the htseq-count command would be run within docker and the redirection would be handled by the host, which has the unfortunate consequence that stderr output would be included in the file. This means that my_counts.txt would also contain diagnostic messages like "1,000,000 SAM lines processed". We work around this by running. [Here](http://www-huber.embl.de/users/anders/HTSeq/doc/count.html) is the htseq-count documentation to read about the options:

```
docker run -i -t --rm -v $(pwd):/data bschiffthaler/ngs bash -c 'htseq-count -f bam -r \
pos -s no -t exon -i Parent /data/202_subset_1_sorted_trimmed_STARAligned.sortedByCoord.out.bam \
/data/Ptrichocarpa_v3.0_210_synthetic-gene-models.gff3 > \
/data/202_subset_1_sorted_trimmed_STAR_HTSeq.txt'
```

## Downstream analysis using RStudio

Finally, the data is ready to be analysed by your differential expression testing tool of choice (e.g.: edgeR, DESeq2, limma, ...). We can therefore start RStudio like this:

```
docker run -p 8787:8787 --rm -v $(pwd):/data bschiffthaler/ngs
```

This spawns an RStudio server, which you can reach by navigating to the host's IP on port 8787 in your web browser i.e.: if my IP is 192.168.1.10 (you can check this by running `ifconfig` on UNIX based systems and `ipconfig` on Windows) I would open my broweser and type:

```
192.168.1.10:8787
```

in the address bar. You can then login with the credentials rstudio as both username and password. On windows and Mac OS, docker runs inside a virtual machine, you therefore need to provide the virtual machine's IP (which you see when you run `boot2docker` or by checking the $DOCKER_HOST environment variable). From within RStudio, you could run:

```
setwd("/data")
dir()
```

You should see your generated data, which you can now further process inside R. **CAVEAT: Save any results etc. only into the folder you have mounted inside the container (in this example "/data") as docker is ephemeral (when run with the '--rm' option) and will not keep any files.**

For the actual analysis, I would recommend the documentation of [DESeq2](http://www.bioconductor.org/packages/release/bioc/html/DESeq2.html), [edgeR](http://master.bioconductor.org/packages/release/bioc/html/edgeR.html) or [limma](http://master.bioconductor.org/packages/release/bioc/html/limma.html). An R transcript of how the example data was analysed is available [here](https://microasp.upsc.se/root/upscb-public/blob/master/projects/Robinson2014/src/Robinson2014-rnaseq-analysis.R).


**References**
----------------------

[1] Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics (Oxford, England), 1–7. doi:10.1093/bioinformatics/btu170

[2] Kopylova, E., Noé, L., & Touzet, H. (2012). SortMeRNA: Fast and accurate filtering of ribosomal RNAs in metatranscriptomic data. Bioinformatics, 28, 3211–3217. doi:10.1093/bioinformatics/bts611

[3] Dobin, A., Davis, C. A., Schlesinger, F., Drenkow, J., Zaleski, C., Jha, S., … Gingeras, T. R. (2013). STAR: Ultrafast universal RNA-seq aligner. Bioinformatics, 29, 15–21. doi:10.1093/bioinformatics/bts635

[4] Li, H., & Durbin, R. (2009). Fast and accurate short read alignment with Burrows-Wheeler transform. Bioinformatics, 25, 1754–1760. doi:10.1093/bioinformatics/btp324

[5] Anders, S., Pyl, P. T., & Huber, W. (2014). HTSeq A Python framework to work with high-throughput sequencing data. bioRxiv. doi:10.1101/002824

[6] Gentleman, R. C., Carey, V. J., Bates, D. M., Bolstad, B., Dettling, M., Dudoit, S., … Zhang, J. (2004). Bioconductor: open software development for computational biology and bioinformatics. Genome Biology, 5, R80. doi:10.1186/gb-2004-5-10-r80

[7] R Core Team. (2014). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria, URL http://www.R–project.org/.

[8] Tuskan, G. A., Difazio, S., Jansson, S., Bohlmann, J., Grigoriev, I., Hellsten, U., … Rokhsar, D. (2006). The genome of black cottonwood, Populus trichocarpa (Torr. & Gray). Science (New York, N.Y.), 313, 1596–1604. doi:10.1126/science.1128691

[9] Robinson, K. M., Delhomme, N., Mähler, N., Schiffthaler, B., Onskog, J., Albrectsen, B. R., … Street, N. R. (2014). Populus tremula (European aspen) shows no evidence of sexual dimorphism. BMC Plant Biology, 14, 276. doi:10.1186/s12870-014-0276-5

[10] Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., … Durbin, R. (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics, 25, 2078–2079. doi:10.1093/bioinformatics/btp352

[11] Li, H. (2011). A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. Bioinformatics, 27, 2987–2993. doi:10.1093/bioinformatics/btr509
