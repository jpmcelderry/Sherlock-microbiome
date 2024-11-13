#!/bin/bash

ml samtools
ml kraken
ml bracken
ml trimmomatic
ml hisat
ml bwa

name=$2
Ssource=$1
db=/lscratch/$SLURM_JOB_ID/k2_db

if [[ $Ssource == "RNA" ]]
then
	input=$3
	confidences=(0.05 0.1 0.15 0.25)
	#db=/data/Sherlock_Lung/JohnMce/k2_pluspf_20230605
	#db=/lscratch/$SLURM_JOB_ID/k2_standardPF_plusTranscriptome
	species=5
	genus=10
	family=10
elif [[ $Ssource == "WGS" ]]
then
	input=$3
	confidences=(0.05 0.1 0.15 0.25)
	#db=/lscratch/$SLURM_JOB_ID/JohnMce/k2_pluspf_20230605
	#db=/data/Sherlock_Lung/JohnMce/k2_standardPF_plusTranscriptome
	species=2
	genus=5
	family=5
elif [[ $Ssource == "16s" ]]
then
	file1=$3
	file2=$4
	confidences=(0.02 0.03 0.04)
	#db=/data/Sherlock_Lung/JohnMce/ncbi_16s
	species=2
	genus=5
	family=5
else
 	echo "unrecognized sequencing type"
  	exit 1
fi

echo ${name}

if [ ! -d ./$name ]
then
	mkdir $name
fi

if [ ! -f ${name}/trimmed1.fq.gz ]
then
	if [[ $Ssource == "RNA" ]]
	then
		##extract both unaligned
		samtools view -@ $SLURM_CPUS_PER_TASK -f 4,8 -F 512,1024 --output-fmt bam ${input}| samtools sort -@ $SLURM_CPUS_PER_TASK -n | samtools fastq -@ $SLURM_CPUS_PER_TASK -s ${name}/${name}_single.fq -1 ${name}/${name}_1.fq -2 ${name}/${name}_2.fq


		##realign with hisat2 and extract both unaligned
		hisat2 -p $SLURM_CPUS_PER_TASK -x /data/Sherlock_Lung/JohnMce/CHM13/chm13v2.0 -1 ${name}/${name}_1.fq -2 ${name}/${name}_2.fq --sensitive | samtools view -@ $SLURM_CPUS_PER_TASK -f 4,8 -F 512,1024 --output-fmt bam | samtools sort -@ $SLURM_CPUS_PER_TASK -n | samtools fastq -@ $SLURM_CPUS_PER_TASK -0 ${name}/cleaned_single1_v2.fq -s ${name}/cleaned_single2_v2.fq -1 ${name}/chm13_unalign1.fq -2 ${name}/chm13_unalign2.fq

	elif [[ $Ssource == "WGS" ]]
	then
		##extract both unaligned
		samtools view -@ $SLURM_CPUS_PER_TASK -f 4,8 -F 512,1024 --output-fmt bam ${input}| samtools sort -@ $SLURM_CPUS_PER_TASK -n | samtools fastq -@ $SLURM_CPUS_PER_TASK -s ${name}/${name}_single.fq -1 ${name}/${name}_1.fq -2 ${name}/${name}_2.fq


		##realign with bwa-mem and extract both unaligned
		bwa mem -t $SLURM_CPUS_PER_TASK /data/Sherlock_Lung/JohnMce/CHM13/chm13v2.0.fa ${name}/${name}_1.fq ${name}/${name}_2.fq | samtools view -@ $SLURM_CPUS_PER_TASK -f 4,8 -F 512,1024 --output-fmt bam | samtools sort -@ $SLURM_CPUS_PER_TASK -n | samtools fastq -@ $SLURM_CPUS_PER_TASK -s ${name}/cleaned_single2_v2.fq -1 ${name}/chm13_unalign1.fq -2 ${name}/chm13_unalign2.fq

	else
		bwa mem -t $SLURM_CPUS_PER_TASK /data/Sherlock_Lung/JohnMce/CHM13/chm13v2.0.fa ${file1} ${file2} | samtools view -@ $SLURM_CPUS_PER_TASK -f 4,8 -F 512,1024 --output-fmt bam | samtools sort -@ $SLURM_CPUS_PER_TASK -n | samtools fastq -@ $SLURM_CPUS_PER_TASK -0 ${name}/cleaned_single1_v2.fq -s ${name}/cleaned_single2_v2.fq -1 ${name}/chm13_unalign1.fq -2 ${name}/chm13_unalign2.fq

	fi

	java -jar $TRIMMOJAR PE -phred33 -trimlog ${name}/trimlog.log ${name}/chm13_unalign1.fq ${name}/chm13_unalign2.fq ${name}/trimmed1.fq ${name}/trimmed_UP1.fq ${name}/trimmed2.fq ${name}/trimmed_UP2.fq \
			ILLUMINACLIP:/usr/local/apps/trimmomatic/0.39/adapters/TruSeq3-PE.fa:2:30:10:2:True SLIDINGWINDOW:4:10 MINLEN:45
	
 	cat ${name}/trimmed_UP1.fq ${name}/trimmed_UP2.fq > ${name}/trimmed_UP.fq

	gzip ${name}/*.fq
fi

for conf in "${confidences[@]}"
do
  kraken2 \
  --threads $SLURM_CPUS_PER_TASK -db ${db} \
  --paired --gzip-compressed --confidence ${conf} \
  --output ${name}/${name}-kraken-conf${conf}-pair.txt \
  --report ${name}/${name}-kraken-report-conf${conf}-pair.txt \
  --minimum-hit-groups 2 \
  ${name}/trimmed1.fq.gz ${name}/trimmed2.fq.gz

  kraken2 \
  --threads $SLURM_CPUS_PER_TASK -db ${db} \
  --gzip-compressed --confidence ${conf} \
  --output ${name}/${name}-kraken-conf${conf}-unpair.txt \
  --minimum-hit-groups 2 \
  ${name}/trimmed_UP.fq.gz

  python ~/KrakenTools-master/make_kreport.py -i <(cat ${name}/${name}-kraken-conf${conf}-pair.txt ${name}/${name}-kraken-conf${conf}-unpair.txt) -t ${db}/ktaxonomy.tsv -o ${name}/${name}-kraken-report-conf${conf}.txt
	# run bracken with shotgun seq parameters
	if [[ $Ssource == "RNA" || $Ssource == "WGS" ]]
	then
    python /data/Sherlock_Lung/JohnMce/est_abundance.py \
    -k ${db}/database100mers.kmer_distrib \
    -i ${name}/${name}-kraken-report-conf${conf}.txt \
    -o ${name}/${name}-braken-conf${conf}.txt \
    --out-report ${name}/${name}-braken-report-conf${conf}.txt \
    -t ${species} -l S
      python /data/Sherlock_Lung/JohnMce/est_abundance.py \
    -k ${db}/database100mers.kmer_distrib \
    -i ${name}/${name}-kraken-report-conf${conf}.txt \
    -o ${name}/${name}-braken-conf${conf}-Genus.txt \
    --out-report ${name}/${name}-braken-report-conf${conf}-Genus.txt \
    -t ${genus} -l G
		python /data/Sherlock_Lung/JohnMce/est_abundance.py \
	 -k ${db}/database100mers.kmer_distrib \
	 -i ${name}/${name}-kraken-report-conf${conf}.txt \
	 -o ${name}/${name}-braken-conf${conf}-Family.txt \
	 --out-report ${name}/${name}-braken-report-conf${conf}-Family.txt \
	 -t ${family} -l F
	 #run bracken with 16s parameters
	else
		python /data/Sherlock_Lung/JohnMce/est_abundance.py \
		-k ${db}/database100mers.kmer_distrib \
		-i ${name}/${name}-kraken-report-conf${conf}.txt \
		-o ${name}/${name}-braken-conf${conf}.txt \
		--out-report ${name}/${name}-braken-report-conf${conf}.txt \
		-t ${species} -l G
     python /data/Sherlock_Lung/JohnMce/est_abundance.py \
    -k ${db}/database100mers.kmer_distrib \
    -i ${name}/${name}-kraken-report-conf${conf}.txt \
    -o ${name}/${name}-braken-conf${conf}-Family.txt \
    --out-report ${name}/${name}-braken-report-conf${conf}-Family.txt \
    -t ${family} -l F
	fi
done

module list 2> ${name}/modules.log
