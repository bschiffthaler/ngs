#!/bin/bash
JBROWSEDIR="/var/www/html/JBrowse-1.11.6/data"
if [ $# -gt 0 ]; then
	INDIR=$@
else
	INDIR="/home/training"
fi
>$JBROWSEDIR/tracks.conf

find $INDIR -name "*.bam" -type f | while read l; do
	chmod o+rw $l
	ln -sf $l ${JBROWSEDIR}/bam
	ln -sf ${l}.bai ${JBROWSEDIR}/bam
	BNAM=$(basename $l)

	echo "[ tracks . ${BNAM//./_}_Alignment ]" >> $JBROWSEDIR/tracks.conf
	echo "storeClass  = JBrowse/Store/SeqFeature/BAM" >> $JBROWSEDIR/tracks.conf
	echo "type        = JBrowse/View/Track/Alignments2" >> $JBROWSEDIR/tracks.conf
	echo "key         = ${BNAM//./_}_Alignment" >> $JBROWSEDIR/tracks.conf
	echo "urlTemplate = bam/$BNAM" >> $JBROWSEDIR/tracks.conf
	echo >> $JBROWSEDIR/tracks.conf
	echo "[ tracks . ${BNAM//./_}_SNP_Histogram ]" >> $JBROWSEDIR/tracks.conf
	echo "storeClass  = JBrowse/Store/SeqFeature/BAM" >> $JBROWSEDIR/tracks.conf
	echo "type        = JBrowse/View/Track/SNPCoverage" >> $JBROWSEDIR/tracks.conf
	echo "key         = ${BNAM//./_}_SNP_Histogram" >> $JBROWSEDIR/tracks.conf
	echo "urlTemplate = bam/$BNAM" >> $JBROWSEDIR/tracks.conf
	echo >> $JBROWSEDIR/tracks.conf
done
