# News

- Nov 15th 2019: Updated software to most recent stable versions (at the time of this update)
- Nov 15th 2019: Changed GateOne (unmaintained) to shellinabox

- I recorded myself doing basic QC on the training data. Check [here](https://www.youtube.com/watch?v=1rNEkWSxB5s) for the video.

# General information

Ready-to-work docker for next generation sequence analysis including binaries:

- Sequence data QC ([FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [MultiQC](http://multiqc.info/)) [15]
- Trimming [(Trimmomatic)](http://www.usadellab.org/cms/?page=trimmomatic) [1]
- rRNA filtering [(SortMeRNA)](http://bioinfo.lifl.fr/RNA/sortmerna/) [2]
- Genome mapping ([STAR](https://github.com/alexdobin/STAR) [3] , [BWA](http://bio-bwa.sourceforge.net/) [4] , [kallisto](https://pachterlab.github.io/kallisto/) [14], [salmon](https://combine-lab.github.io/salmon/)[16])
- Feature Summarisation [(HTSeq)](http://www-huber.embl.de/HTSeq/doc/overview.html) [5]
- File manipulation and exploration [(samtools,htslib,bcftools)](http://www.htslib.org/) [8],[9]
- Alignment visualisation ([JBrowse](http://jbrowse.org/)) [10]
- Peak calling ([MACS2](http://liulab.dfci.harvard.edu/MACS/))  [11]
- Sequence data analysis ([Useq](http://useq.sourceforge.net/)) [12]
- Binding site determination ([SISSRs](http://www.rajajothi.com/sissrs/)) [13]

For downstream analysis, this docker is based on bioconductor/release_core2 [6] , which contains all the most commonly used downstream analysis tools implemented in R [7] .

The ":with-data" tagged image used to contain a set of training data which is commonly used in our RNA-Seq training courses. Due to size concerns of the docker image, this tag is no longer available.

The source (+Dockerfile) which was used to build this container, is [here](https://github.com/bschiffthaler/ngs) on GitHub!

# Common use cases

Even if meant for training purposes, this Docker should allow a scientist to carry out an RNA-Seq based differential expression study as described [in this protocol](http://www.epigenesys.eu/images/stories/protocols/pdf/20150303161357_p67.pdf).

For reasons of brevity, I will not be going into great detail (see the protocol), but I will give a short overview of common use cases in differential expression analysis.

## Index

1. Basics
   - Mounting host directories inside the docker container
2. Technical quality control with FastQC
3. rRNA sorting with SortMeRNA
4. Quality trimming with Trimmomatic
5. Mapping to the reference genome with STAR
   - Creating a STAR genome
   - Mapping my reads
6. Visualisation of the alignments in JBrowse
7. Feature summarisation using HTSeq
8. Downstream analysis using RStudio
9. References

## Basics

### Mounting host directories inside the docker container

In order to analyse your data from within the docker, you need to mount the directory containing it in the docker container. In this example, I have a FASTQ file called my_sample.fq.gz in my current working directory and would like to analyse it with software from my docker. I therefore mount my current working directory in the `/data` folder within the container. The `--rm` flag tells docker to clean up the process after it's done.:

```
docker run --rm -v $(pwd):/data bschiffthaler/ngs zcat /data/my_sample.fq.gz | head
```

The argument `-v $(pwd):/data` tells docker to mount the current working directory (from the `pwd` command) in the folder called `/data` within the container. Then, the `zcat` command is executed, which decompresses the file and streams the contents to the terminal. We pipe (using `|`) the output to the `head` command so that we only read the first few lines. Of course, this example is almost entirely non practical as you could do all of this without having the docker image.

## Technical quality control with FastQC

As a small prologue to this section I would like to differentiate between technical and biological QA. As an example, technical QA would show sequencing errors, PCR over-amplification and similar issues, while a sample swap during the RNA extraction would not show up during this part. For the lack of a better term, we call methods to determine the latter type of issues "biological QA". This part is about technical QA.

After you receive your sequencing data, and after every step that modifies it (trimming, rRNA sorting) QA should be done. Here, we will be using Simon Andrews' [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/):

```
#Here I assume that I have a FASTQ file called 202_subset_1.fq.gz in my current working directory
docker run --rm -v $(pwd):/data bschiffthaler/ngs fastqc /data/202_subset_1.fq.gz
```

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

## Alignment visualisation using JBrowse

Unlike the previous commands, we need to keep the docker container alive and interactive, as JBrowse requires several steps to be functional. We will further need to map port 80 from the docker to the host in order to view alignments in apache. The option `-p 80:80` maps the port from the container to the host. Once we start apache, you can then open a web browser and navigate to your IP (on Linux hosts you can simply type "localhost" in the browser's address bar, on OSX and Windows hosts, you can find your docker's IP address if you run `boot2docker ip` from the boot2docker command. Look [here](https://docs.docker.com/userguide/usingdocker/#viewing-our-web-application-container) for more info.) The other options `-ti` tell docker that we want an interactive session, and that we want a (pseudo-)terminal to enter commands.

```
docker run -p 80:80 -v $(pwd):/data -ti bschiffthaler/ngs bash -l
```

First, we start by adding our reference(s) to JBrowse. This is done with perl scripts which are delivered with JBrowse. The reference sequence in FASTA format and the reference annotation in GFF3 format are in my current working directory, as are the alignment results in BAM format.

Let's start by adding the FASTA and GFF3 files as tracks in JBrowse:

```
cd /var/www/html/JBrowse-1.11.6
bin/prepare-refseqs.pl --fasta /data/Ptrichocarpa_v3.0_210.fa
bin/flatfile-to-json.pl --gff /data/Ptrichocarpa_210_gene_exons.gff3 --trackType CanvasFeatures --trackLabel P.trichocarpa_v3.0_gene_exons
```

This has now added the reference sequence and annotation as tracks in JBrowse. Now let's add our alignments. I added a script to this docker, which automatically adds a read alignment viewer track as well as a SNP and coverage histogram track to the file `/var/www/html/JBrowse-1.11.6/data/tracks.conf`. My alignment files are in a sub-folder within my current working directory. Before you run the script to add tracks, please make sure that your alignment file is BAM formatted, sorted and indexed. If you ran `STAR` as described above, you should already have a sorted BAM file. If you haven't you can sort the file:

```
samtools sort /data/<my_BAM_file.bam> -o /data/<my_sorted_BAM_file.bam>
```

Then index the file

```
samtools index /data/<my_sorted_BAM_file.bam>
```

After that is done, we can add the tracks to JBrowse

```
add_JBrowse_tracks.sh /data
```

Finally, we can start the web server apache:

```
service apache2 start
```

Now, if we open a browser on the host and navigate to the aforementioned IP address (or localhost), we should see a link to JBrowse as well as the training user's home folder. Click on the JBrowse link, which will start JBrowse and browse your alignments!

I added a screenshot [here](https://www.dropbox.com/s/mw2uvugrh06y44y/Screenshot%202015-04-10%2002.43.22.png?dl=0) of what you should see!

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

## **References**

[1] Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics (Oxford, England), 1–7. doi:10.1093/bioinformatics/btu170

[2] Kopylova, E., Noé, L., & Touzet, H. (2012). SortMeRNA: Fast and accurate filtering of ribosomal RNAs in metatranscriptomic data. Bioinformatics, 28, 3211–3217. doi:10.1093/bioinformatics/bts611

[3] Dobin, A., Davis, C. A., Schlesinger, F., Drenkow, J., Zaleski, C., Jha, S., … Gingeras, T. R. (2013). STAR: Ultrafast universal RNA-seq aligner. Bioinformatics, 29, 15–21. doi:10.1093/bioinformatics/bts635

[4] Li, H., & Durbin, R. (2009). Fast and accurate short read alignment with Burrows-Wheeler transform. Bioinformatics, 25, 1754–1760. doi:10.1093/bioinformatics/btp324

[5] Anders, S., Pyl, P. T., & Huber, W. (2014). HTSeq A Python framework to work with high-throughput sequencing data. bioRxiv. doi:10.1101/002824

[6] Gentleman, R. C., Carey, V. J., Bates, D. M., Bolstad, B., Dettling, M., Dudoit, S., … Zhang, J. (2004). Bioconductor: open software development for computational biology and bioinformatics. Genome Biology, 5, R80. doi:10.1186/gb-2004-5-10-r80

[7] R Core Team. (2014). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria, URL http://www.R–project.org/.

[8] Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., … Durbin, R. (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics, 25, 2078–2079. doi:10.1093/bioinformatics/btp352

[9] Li, H. (2011). A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. Bioinformatics, 27, 2987–2993. doi:10.1093/bioinformatics/btr509

[10] Skinner, M. E., Uzilov, A. V., Stein, L. D., Mungall, C. J., & Holmes, I. H. (2009). JBrowse: A next-generation genome browser. Genome Research, 19, 1630–1638. doi:10.1101/gr.094607.109

[11] Zhang, Y., Liu, T., Meyer, C. A., Eeckhoute, J., Johnson, D. S., Bernstein, B. E., … Liu, X. S. (2008). Model-based analysis of ChIP-Seq (MACS). Genome Biology, 9(9), R137. doi:10.1186/gb-2008-9-9-r137

[12] Nix, D. A., Courdy, S. J., & Boucher, K. M. (2008). Empirical methods for controlling false positives and estimating confidence in ChIP-Seq peaks. BMC Bioinformatics, 9, 523. doi:10.1186/1471-2105-9-523

[13] Jothi, R., Cuddapah, S., Barski, A., Cui, K., & Zhao, K. (2008). Genome-wide identification of in vivo protein-DNA binding sites from ChIP-Seq data. Nucleic Acids Research, 36(16), 5221–5231. doi:10.1093/nar/gkn488

[14] Bray, N. L., Pimentel, H., Melsted, P., & Pachter, L. (2015). Near-optimal RNA-Seq quantification. aRxiv. http://doi.org/arXiv:1505.02710

[15] Ewels, P., Magnusson, M., Lundin, S., & K??ller, M. (2016). MultiQC: Summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 32(19), 3047–3048. http://doi.org/10.1093/bioinformatics/btw354

[16] Patro, R., Duggal, G., Love, M. I., Irizarry, R. A., & Kingsford, C. (2017). Salmon provides fast and bias-aware quantification of transcript expression. Nature Methods.
