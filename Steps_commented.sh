#!/bin/bash

# This is not the script itself, but a summary of the main command lines used to make the presentation easier.
# contact Aurelie PETICCA if you have any questions 
#
# Explanation of variables used  
# $adapters_file = file containing the sequences of adapters used for sequencing  
# $blacklisted_file = file containing black listed sequences for peak calling (see https://doi-org.insb.bib.cnrs.fr/10.1038/s41598-019-45839-z)
# $folder_out = name of writing folder / output folder
# $genome_file = file containing the genome sequences (fasta file)
# $gtf_file = file containing the annotation fo the genome (gtf file)
# $name_output = name of writing file / output file
# $number_reads_to_downsample = maximum number of reads for each replicate of a sample (downsampling), calculate using a homemade script 
# $peaks_file = MACS output file containing the peaks called (bed file)
# $reads1 / $reads2 = files containing raw sequencing data (here Illumina R1 and R2) 
# $reads1_trim / $reads_trim = files containing trimmed sequencing data (here Illumina R1 and R2) 
# $region_to_check = files containing the regions where ChIP values are to be compared (bed format)
# $replicate_name = name of the replicate of a sample studied
# $sample_name = name of the sample studied ($sample_name.bam = merged of downsampled bam replicates) 
# $SLURM_CPUS_PER_TASK = Number of CPUs allocated to each task (SLURM Workload Manager)  

 
#Cutadapt 
cutadapt \
    -j $SLURM_CPUS_PER_TASK \
    -o $folder_out \
    -p $replicate_name \
    -a file:$adapters_file \
    -A file:$adapters_file \
    --minimum-length 75:75 \
    -q 20 \
    --pair-filter=any \
    --quality-base=33 \
    $reads1 $reads2 
    # -a file:$adapters_file : adaptaters 5' for R1
    # -A file:$adapters_file : adaptaters 5' for R2
    # --minimum-length 75:75 : minimum length for R1:R2
    # -q 20 : minimum quality phred
    # --pair-filter=any : discard if one read of the pair isn't ok
    # --quality-base=33 : output in phred33

#Bowtie2
bowtie2 \
    -p $SLURM_CPUS_PER_TASK  \
    --very-sensitive \
    --phred33 \
    --no-mixed \
    --no-discordant \
    --dovetail \
    -x $genome_file.index \
    -1 $reads1_trim -2 $reads2_trim \
    -S $replicate_name.sam
    # --very-sensitive : more sensitive and more accurate, same as running with options: -D 20 -R 3 -N 0 -L 20 -i S,1,0.50 
    # --phred33 : to be sure that phred score taken is 33
    # --no-unal : don't output in SAM the unaligned PE , NOT USED HERE FOR STATISTICS PURPOSE AFTER
    # --no-mixed : only consider alignment status of pairs per se
    # --no-discordant : the readings are directed towards each other
    # --dovetail : the reads can overloap themself 

#Samtools view/sort
samtools view \
    -F 1804 \
    -f 2 \
    -q 20 \
    -bS $replicate_name.sam \
	> $replicate_name.bam 
    # -q 20 : have mapping quality >= 20 , https://doi-org.insb.bib.cnrs.fr/10.1038/s41467-022-35384-1
    # -F 1804 : remove reads unmapped, mate unmapped, not primary alignment, reads failing platform, duplicates
    # -f 2 : retain properly paired reads 
samtools sort -O BAM \
    -o $replicate_name.sort.bam \
	$replicate_name.bam

#Picard MarkDuplicate
picard MarkDuplicates \
    --INPUT $replicate_name.sort.bam \
    --OUTPUT $replicate_name.bam \
    --REMOVE_DUPLICATES \
    --METRICS_FILE $replicate_name.picstats 
    # --REMOVE_DUPLICATES : delete duplicates in the output file, mandatory for ChIP-Seq analyses 

#Bedtools intersect 
bedtools intersect \
    -abam $replicate_name.sort.bam \
    -b $blacklisted_file \
    -v \
    -sorted \
    > $replicate_name.bam
samtools sort \
    -O BAM \
    -o $replicate_name.sort.bam \
    $replicate_name.bam

#Samtools downsampling, merge and sort
fraction_to_keep=$(samtools idxstats $replicate_name.sort.bam \
    | cut -f3 | awk -v ct=$number_reads_to_downsample 'BEGIN {total=0} {total += $1} END {print ct/total}')
samtools view \
    -b \
    --subsample $fraction_to_keep \
    -@ $SLURM_CPUS_PER_TASK \
    > $sample_name.bam
samtools sort \
    -O BAM \
    -o $sample_name.sort.bam \
    $sample_name.bam

#MACS+Homer
macs3 callpeak \
    -t $sample_name.sort.bam \
    --nolambda \
    -f BAMPE \
    -g 2495461690 \
    -n $sample_name \
    -B \
    -q 0.01 \
    --nomodel \
    --shift 0 \
    --extsize 147 \
    --call-summits \
    --outdir $folder_out
#OR
macs3 callpeak \
    -t $sample_name.sort.bam \
    -c $control_sample_name.sort.bam \
    -f BAMPE \
    -g 2495461690 \
    -n $sample_name \
    -B \
    -q 0.01 \
    --nomodel \
    --shift 0 \
    --extsize 147 \
    --call-summits \
    --outdir $folder_out
    # -g found here : https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html
    # -q 0.01 for IF OR 0.05 for Marked histo, both were used in this study
    # --shift 0 : because BAMPE (see MACS manual : https://macs3-project.github.io/MACS/docs/callpeak.html)
annotatePeaks.pl \
    $peaks_file \
    $genome_file \
    -gtf $gtf_file \
    > $sample_name.narrowPeaks.annotated.bed

#Deeptools normalization
#Before
bamCompare \
    --bamfile1 $sample_name.sort.bam \
    --bamfile2 $control_sample_name.sort.bam \
    --outFileFormat bigwig \
    --outFileName $sample_name.log2ratio.bw \
    --effectiveGenomeSize 2495461690 \
    --normalizeUsing RPGC \
    --centerReads \
    --extendReads \
    --smoothLength 3 \
    --binSize 10
    # --effectiveGenomeSize found here : https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html
    # --binSize 10 : smaller bin size for higher resolution
    # --normalizeUsing RPGC : reads per genomic content (1x normalization) 
    # --centerReads : read centered at fragment length defined by the two ends of the fragment, sharper signal around enriched regions
    # --extendReads : reads with mates are always extended to match the fragment size defined by the two read mates
    # Have also tried : --scaleFactorsMethod None and --normalizeUsing BPM
#Now
bamCoverage \
    -b $sample_name.sort.bam \
    --outFileFormat bigwig \
    --outFileName $sample_name.rpgc.bw \
    --normalizeUsing RPGC \
    --effectiveGenomeSize 2495461690 \
    --MNase \
    --numberOfProcessors $SLURM_CPUS_PER_TASK \
    --minFragmentLength 80 \
    --maxFragmentLength 200 \
    --smoothLength 80 
    # --effectiveGenomeSize found here : https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html
    # --normalizeUsing RPGC : RPGC (per bin) = number of reads per bin / scaling factor for 1x average coverage
    # --maxFragmentLength 200 : distribution of the fragments, 167 will remove ~50% of the data 
    # --smoothLength 80 : test running
    # --MNase : Because we work with an MNase dataset, it helps position the nucleosome. 

#Deeptools plots
computeMatrix reference-point \
    -S $sample_name.rpgc.bw \
    -R $region_to_check \
    --beforeRegionStartLength 1000 \
    --afterRegionStartLength 1000 \
    --referencePoint center \
    --sortRegions keep \
    --numberOfProcessors $SLURM_CPUS_PER_TASK \
    --smartLabels \
    --metagene \
    --outFileName $name_output.$region_to_check.bw.gz
    # --referencePoint center or TSS : depends the region to check
plotProfile -m $name_output.$region_to_check.bw.gz \
    -o $name_output.$region_to_check.profile.svg \
    --plotFileFormat svg \
    --perGroup \
    --plotType lines \
    --numPlotsPerRow 1 \
    --plotHeight 9 \
    --plotWidth 7 \
    --samplesLabel Label1 Label2 Label3 ... \
    --legendLocation best \
    --plotTitle "Title"

############################################################################
# Supplementary data, log2ratio example
bigwigCompare -b1 HA-CT20-WT.sort.bam.rpgc.bw \ 
  -b2 H3-12-CT20-WT.sort.bam.rpgc.bw \ 
  --operation log2 \ 
  --pseudocount 1 \ 
  -o log2_H3.3_vs_H3.12_CT20.bw

computeMatrix reference-point \
  -S log2_H3.3_vs_H3.12_CT8.bw log2_H3.3_vs_H3.12_CT20.bw \
  -R Mus_musculus.GRCm39.112.genes.bed \
  --beforeRegionStartLength 1000 \
  --afterRegionStartLength 1000 \
  --referencePoint TSS \
  --smartLabels --metagene \
  -o log2_H3.3_H3.12_CT8_vs_CT20.gz
  
plotProfile \
  -m log2_H3.3_H3.12_CT8_vs_CT20.gz \
  --perGroup \
  --plotType lines \
  --plotFileFormat svg \
  --legendLocation best \
  --outFileName log2_H3.3_H3.12_CT8_vs_CT20.profile.svg \
  --outFileNameData log2_H3.3_H3.12_CT8_vs_CT20.profile.tab \
  --colors blue red --plotHeight 7 --plotWidth 9 \
  --samplesLabel "log2(H3.3/H3.1-2) CT8" "log2(H3.3/H3.1-2) CT20" \
  --yMin 0.02 --yMax 0.27
