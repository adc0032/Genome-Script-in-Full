# Genome-Script-in-Full

##Indexing the Reference: Step One 

```ruby
#! /bin/bash

source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load bwa/0.7.12
#module load samtools/1.2
#module load picard/1.79

#index sacCer3 genome file
bwa index -p sacCer3 ~/Group_Project/sacCer3.masked.fa

#samtools faidx ./sacCer3/sacCer3.fa.masked
#gunzip ~/class_shared/sacCer3/sacCer3.fa.masked.gz
#mv ~/class_shared/sacCer3/sacCer3.fa.masked ~/class_shared/sacCer3/sacCer3.masked.fa

#samtools faidx ~/class_shared/sacCer3/sacCer3.masked.fa
#cd ~/class_shared/sacCer3/
#java -Xms2g -Xmx4g -jar /opt/asn/apps/picard_1.79/picard-tools-1.79/CreateSequenceDictionary.jar REFERENCE=sacCer3.masked.fa OUTPUT=sacCer3.masked.dict
```

##Download the Sample Files: Step Two 
```ruby
#!/bin/sh

source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load sra/2.3.5

fastq-dump --origfmt --gzip -I SRR1693723
fastq-dump --origfmt --gzip -I SRR1693724
fastq-dump --origfmt --gzip -I SRR1693728


```
##Split Lanes: Step Three
```ruby
#!/bin/bash
#copy your fastq file into scratch
cp ~/Group_Project/Part2/Original_downloads/SRR1693723.fastq.gz /scratch/aubcls35/
cd /scratch/aubcls35/
gunzip SRR1693723.fastq.gz
awk 'BEGIN {FS = ":"} {lane=$3.$4 ; print > "SRR1693723."lane".fastq" ; for (i = 1; i <= 3; i++) {getline ; print > "SRR1693723."lane".fastq"}}' SRR1693723.fastq
gzip SRR1693723.*.fastq
cp SRR1693723.*.fastq.gz ~/Group_Project
rm SRR1693723.fastq

cp ~/Group_Project/Part2/Original_downloads/SRR1693724.fastq.gz /scratch/aubcls35/
cd /scratch/aubcls35/
gunzip SRR1693724.fastq.gz
awk 'BEGIN {FS = ":"} {lane=$3.$4 ; print > "SRR1693724."lane".fastq" ; for (i = 1; i <= 3; i++) {getline ; print > "SRR1693724."lane".fastq"}}' SRR1693724.fastq
gzip SRR1693724.*.fastq
cp SRR1693724.*.fastq.gz ~/Group_Project
rm SRR1693724.fastq

cp ~/Group_Project/Part2/Original_downloads/SRR1693728.fastq.gz /scratch/aubcls35/
cd /scratch/aubcls35/
gunzip SRR1693728.fastq.gz
awk 'BEGIN {FS = ":"} {lane=$3.$4 ; print > "SRR1693728."lane".fastq" ; for (i = 1; i <= 3; i++) {getline ; print > "SRR1693728."lane".fastq"}}' SRR1693728.fastq
gzip SRR1693728.*.fastq
cp SRR1693728.*.fastq.gz ~/Group_Project
rm SRR1693728.fastq
```
##BWA Script (Alignment): Step Four
```ruby
#! /bin/bash

#aligns
#M= picard tools compatability
#-v= verbose
#T=Number of threads
#order: reference, input, pipe to output
#sort
#Sbu sam to bam file

module load samtools
module load bwa

#mkdir ~/Group_Project/Bam_files

bwa mem -M  -v 2 -t 6 -R "@RG\tID:Rep9Gen6\tSM:Rep9Gen6Lane1\tPL:illumina\tPU:_1\tLB:SRR1693723_1" sacCer3 /home/aubcls35/Group_Project/SRR1693723.C1630ACXX6.fastq.gz | samtools view -Sb | samtools sort > SRR1693723_1.sorted.bam
bwa mem -M  -v 2 -t 6 -R "@RG\tID:Rep9Gen6\tSM:Rep9Gen6Lane2\tPL:illumina\tPU:_2\tLB:SRR1693723_2" sacCer3 /home/aubcls35/Group_Project/SRR1693723.C1AC9ACXX5.fastq.gz | samtools view -Sb | samtools sort > SRR1693723_2.sorted.bam
bwa mem -M  -v 2 -t 6 -R "@RG\tID:Rep9Gen6\tSM:Rep9Gen6Lane3\tPL:illumina\tPU:_3\tLB:SRR1693723_3" sacCer3 /home/aubcls35/Group_Project/SRR1693723.C1AC9ACXX6.fastq.gz | samtools view -Sb | samtools sort > SRR1693723_3.sorted.bam
bwa mem -M  -v 2 -t 6 -R "@RG\tID:Rep9Gen6\tSM:Rep9Gen6Lane1\tPL:illumina\tPU:_1\tLB:SRR1693724_1" sacCer3 /home/aubcls35/Group_Project/SRR1693724.C1630ACXX6.fastq.gz | samtools view -Sb | samtools sort > SRR1693724_1.sorted.bam
bwa mem -M  -v 2 -t 6 -R "@RG\tID:Rep9Gen6\tSM:Rep9Gen6Lane2\tPL:illumina\tPU:_2\tLB:SRR1693724_2" sacCer3 /home/aubcls35/Group_Project/SRR1693724.C1AC9ACXX5.fastq.gz | samtools view -Sb | samtools sort > SRR1693724_2.sorted.bam
bwa mem -M  -v 2 -t 6 -R "@RG\tID:Rep9Gen6\tSM:Rep9Gen6Lane3\tPL:illumina\tPU:_3\tLB:SRR1693724_3" sacCer3 /home/aubcls35/Group_Project/SRR1693724.C1AC9ACXX6.fastq.gz | samtools view -Sb | samtools sort > SRR1693724_3.sorted.bam
bwa mem -M  -v 2 -t 6 -R "@RG\tID:Rep9Gen6\tSM:Rep9Gen6Lane1\tPL:illumina\tPU:_1\tLB:SRR1693728_1" sacCer3 /home/aubcls35/Group_Project/SRR1693728.C1AC9ACXX2.fastq.gz | samtools view -Sb | samtools sort > SRR1693728_1.sorted.bam
bwa mem -M  -v 2 -t 6 -R "@RG\tID:Rep9Gen6\tSM:Rep9Gen6Lane2\tPL:illumina\tPU:_2\tLB:SRR1693728_2" sacCer3 /home/aubcls35/Group_Project/SRR1693728.C1AC9ACXX3.fastq.gz | samtools view -Sb | samtools sort > SRR1693728_2.sorted.bam
bwa mem -M  -v 2 -t 6 -R "@RG\tID:Rep9Gen6\tSM:Rep9Gen6Lane3\tPL:illumina\tPU:_3\tLB:SRR1693728_3" sacCer3 /home/aubcls35/Group_Project/SRR1693728.D19EFACXX1.fastq.gz | samtools view -Sb | samtools sort > SRR1693728_3.sorted.bam

#index the samples above
samtools index -bc SRR1693723_1.sorted.bam
samtools index -bc SRR1693723_2.sorted.bam
samtools index -bc SRR1693723_3.sorted.bam
samtools index -bc SRR1693724_1.sorted.bam
samtools index -bc SRR1693724_2.sorted.bam
samtools index -bc SRR1693724_3.sorted.bam
samtools index -bc SRR1693728_1.sorted.bam
samtools index -bc SRR1693728_2.sorted.bam
samtools index -bc SRR1693728_3.sorted.bam
```
##Merging: Step Five
```ruby
# !/bin/bash


#Finally merge these together and sort and index the merged bam file to assess the quality using samtools
#have the indexed files in the same directory as the bam files
module load samtools

samtools merge -r SRR1693723.sorted.bam SRR1693723_1.sorted.bam SRR1693723_2.sorted.bam SRR1693723_3.sorted.bam
samtools merge -r SRR1693724.sorted.bam SRR1693724_1.sorted.bam SRR1693724_2.sorted.bam SRR1693724_3.sorted.bam
samtools merge -r SRR1693728.sorted.bam SRR1693728_1.sorted.bam SRR1693728_2.sorted.bam SRR1693728_3.sorted.bam

samtools sort SRR1693723.sorted.bam > SRR1693723.merged.final.bam
samtools sort SRR1693724.sorted.bam > SRR1693724.merged.final.bam
samtools sort SRR1693728.sorted.bam > SRR1693728.merged.final.bam

samtools index -bc SRR1693723.merged.final.bam
samtools index -bc SRR1693724.merged.final.bam
samtools index -bc SRR1693728.merged.final.bam
```

##Statistics: Step Six 

```ruby
module load samtools
samtools flagstat SRR1693723.merged.final.bam
samtools flagstat SRR1693724.merged.final.bam
samtools flagstat SRR1693728.merged.final.bam

samtools depth SRR1693723.merged.final.bam > SRR1693723.depth.out
awk '{sum+=$3} END { print "Average = ",sum/12071326}' SRR1693723.depth.out > SRR723.depth.text
awk '{sum+=$3; sumsq+=$3*$3} END {print "Stdev = ",sqrt(sumsq/12071326 - (sum/12071326)**2)}' SRR1693723.depth.out >> SRR723.depth.text

samtools depth SRR1693724.merged.final.bam > SRR1693724.depth.out
awk '{sum+=$3} END { print "Average = ",sum/12071326}' SRR1693724.depth.out > SRR724.depth.text
awk '{sum+=$3; sumsq+=$3*$3} END {print "Stdev = ",sqrt(sumsq/12071326 - (sum/12071326)**2)}' SRR1693724.depth.out >> SRR724.depth.text

samtools depth SRR1693728.merged.final.bam > SRR1693728.depth.out
awk '{sum+=$3} END { print "Average = ",sum/12071326}' SRR1693728.depth.out > SRR728.depth.text
awk '{sum+=$3; sumsq+=$3*$3} END {print "Stdev = ",sqrt(sumsq/12071326 - (sum/12071326)**2)}' SRR1693728.depth.out >> SRR728.depth.text
```
##Picard and GATK: Step Seven

```javascript
! #/bin/bash

module load picard
module load gatk
module load java

#sort by genomic coordinate (sorted directory)
java -jar /opt/asn/apps/picard_1.79/picard-tools-1.79/SortSam.jar INPUT=SRR1693723.merged.final.bam OUTPUT=SRR1693723.picard.sorted.bam SO=coordinate
java -jar /opt/asn/apps/picard_1.79/picard-tools-1.79/SortSam.jar INPUT=SRR1693724.merged.final.bam OUTPUT=SRR1693724.picard.sorted.bam SO=coordinate
java -jar /opt/asn/apps/picard_1.79/picard-tools-1.79/SortSam.jar INPUT=SRR1693728.merged.final.bam OUTPUT=SRR1693728.picard.sorted.bam SO=coordinate

#mark duplicates (sorted_markdup directory)
java -Xms2g -Xmx10g -jar /opt/asn/apps/picard_1.79/picard-tools-1.79/MarkDuplicates.jar INPUT=SRR1693723.picard.sorted.bam OUTPUT=SRR1693723.picard.sorted.markdup.bam REMOVE_DUPLICATES=false METRICS_FILE=SRR1693723-sort-dup.dup_metrics
java -Xms2g -Xmx10g -jar /opt/asn/apps/picard_1.79/picard-tools-1.79/MarkDuplicates.jar INPUT=SRR1693724.picard.sorted.bam OUTPUT=SRR1693724.picard.sorted.markdup.bam REMOVE_DUPLICATES=false METRICS_FILE=SRR1693724-sort-dup.dup_metrics
java -Xms2g -Xmx10g -jar /opt/asn/apps/picard_1.79/picard-tools-1.79/MarkDuplicates.jar INPUT=SRR1693728.picard.sorted.bam OUTPUT=SRR1693728.picard.sorted.markdup.bam REMOVE_DUPLICATES=false METRICS_FILE=SRR1693728-sort-dup.dup_metrics

#add read groups to file head (only if not done already using samtools)
java -jar /opt/asn/apps/picard_1.79/picard-tools-1.79/AddOrReplaceReadGroups.jar INPUT=SRR1693723.picard.sorted.markdup.bam RGPL=Illumina RGPU=3 RGLB=SRR1693723 RGSM=Rep9Gen6 RGID=SRR1693723 OUTPUT=SRR1693723.sorted.markdup.wHeader.bam
java -jar /opt/asn/apps/picard_1.79/picard-tools-1.79/AddOrReplaceReadGroups.jar INPUT=SRR1693724.picard.sorted.markdup.bam RGPL=Illumina RGPU=4 RGLB=SRR1693724 RGSM=Rep9Gen6 RGID=SRR1693724 OUTPUT=SRR1693724.sorted.markdup.wHeader.bam
java -jar /opt/asn/apps/picard_1.79/picard-tools-1.79/AddOrReplaceReadGroups.jar INPUT=SRR1693728.picard.sorted.markdup.bam RGPL=Illumina RGPU=8 RGLB=SRR1693728 RGSM=Rep9Gen6 RGID=SRR1693728 OUTPUT=SRR1693728.sorted.markdup.wHeader.bam

#make bam index (sorted_markdup directory)
java -Xms2G -Xmx4g -jar /opt/asn/apps/picard_1.79/picard-tools-1.79/BuildBamIndex.jar INPUT=SRR1693723.picard.sorted.markdup.bam OUTPUT=SRR1693723.picard.sorted.markdup.bam.bai
java -Xms2G -Xmx4g -jar /opt/asn/apps/picard_1.79/picard-tools-1.79/BuildBamIndex.jar INPUT=SRR1693724.picard.sorted.markdup.bam OUTPUT=SRR1693724.picard.sorted.markdup.bam.bai
java -Xms2G -Xmx4g -jar /opt/asn/apps/picard_1.79/picard-tools-1.79/BuildBamIndex.jar INPUT=SRR1693728.picard.sorted.markdup.bam OUTPUT=SRR1693728.picard.sorted.markdup.bam.bai

#target regions (intervals directory)
java -Xms2g -Xmx4g -jar /opt/asn/apps/gatk_3.4-46/GenomeAnalysisTK.jar -R sacCer3.masked.fa -T RealignerTargetCreator -I SRR1693723.picard.sorted.markdup.bam -o SRR1693723.forIndelRealigner.intervals
java -Xms2g -Xmx4g -jar /opt/asn/apps/gatk_3.4-46/GenomeAnalysisTK.jar -R sacCer3.masked.fa -T RealignerTargetCreator -I SRR1693724.picard.sorted.markdup.bam -o SRR1693724.forIndelRealigner.intervals
java -Xms2g -Xmx4g -jar /opt/asn/apps/gatk_3.4-46/GenomeAnalysisTK.jar -R sacCer3.masked.fa -T RealignerTargetCreator -I SRR1693728.picard.sorted.markdup.bam -o SRR1693728.forIndelRealigner.intervals

#indel realignment (Final_realigned directory)

java -Xms2g -Xmx10g -jar /opt/asn/apps/gatk_3.4-46/GenomeAnalysisTK.jar -R sacCer3.masked.fa -T IndelRealigner -I SRR1693723.picard.sorted.markdup.bam -o SRR1693723.sorted.markdup.realigned.bam -targetIntervals SRR1693723.forIndelRealigner.intervals
java -Xms2g -Xmx10g -jar /opt/asn/apps/gatk_3.4-46/GenomeAnalysisTK.jar -R sacCer3.masked.fa -T IndelRealigner -I SRR1693724.picard.sorted.markdup.bam -o SRR1693724.sorted.markdup.realigned.bam -targetIntervals SRR1693724.forIndelRealigner.intervals
java -Xms2g -Xmx10g -jar /opt/asn/apps/gatk_3.4-46/GenomeAnalysisTK.jar -R sacCer3.masked.fa -T IndelRealigner -I SRR1693728.picard.sorted.markdup.bam -o SRR1693728.sorted.markdup.realigned.bam -targetIntervals SRR1693728.forIndelRealigner.intervals

#recalibrate bases (we skipped this in our analysis) but Example below
java -jar GenomeAnalysisTK.jar -T PrintReads -R sacCer3.masked.fa -I YeastGenome.realigned.bam -BQSR YeastGen__recalibration_report.grp -o YeastGenome.realigned.bqsr.bam
```











