#!/bin/bash

ml samtools
ml kraken
ml bracken
ml trimmomatic
ml hisat

name=$2
Ssource=$1

if [[ $Ssource == "RNA" ]]
then
	input=$3
	confidences=(0.1 0.15 0.25)
	#db=/data/Sherlock_Lung/JohnMce/k2_pluspf_20230605
	db=/data/Sherlock_Lung/JohnMce/k2_standardPF_plusTranscriptome
	species=5
	genus=10
	family=10
elif [[ $Ssource == "DNA" ]]
then
	input=$3
	confidences=(0.1 0.15 0.25)
	#db=/data/Sherlock_Lung/JohnMce/k2_pluspf_20230605
	db=/data/Sherlock_Lung/JohnMce/k2_standardPF_plusTranscriptome
	species=2
	genus=5
	family=5
elif [[ $Ssource == "16s" ]]
then
	file1=$3
	file2=$4
	confidences=(0.02 0.03 0.04)
	db=/data/Sherlock_Lung/JohnMce/ncbi_16s
	species=2
	genus=5
	family=5
fi

echo ${name}

if [ ! -d ./$name ]
then
	mkdir $name
fi

pigz -p $SLURM_CPUS_PER_TASK ${name}/*.fq
if [ ! -f ${name}/${name}_paired1.fq.gz ]
then
	if [ ! -f ${name}/paired1.fq.gz ]
	then
		if [[ $Ssource == "16s" ]]
		then
		java -jar $TRIMMOJAR PE -phred33 -trimlog ${name}/trimlog.log ${file1} ${file2} ${name}/paired1.fq ${name}/unpaired1.fq ${name}/paired2.fq ${name}/unpaired2.fq \
		ILLUMINACLIP:/usr/local/apps/trimmomatic/0.39/adapters/TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:40 SLIDINGWINDOW:4:10
		else
		samtools fastq -@ $SLURM_CPUS_PER_TASK -f 4 -F 512,1024 -0 /dev/null -s /dev/null -1 ${name}/${name}_1.fastq -2 ${name}/${name}_2.fastq ${input}
		java -jar $TRIMMOJAR PE -phred33 -trimlog ${name}/trimlog.log ${name}/${name}_1.fastq ${name}/${name}_2.fastq ${name}/paired1.fq ${name}/unpaired1.fq ${name}/paired2.fq ${name}/unpaired2.fq \
			ILLUMINACLIP:/usr/local/apps/trimmomatic/0.39/adapters/TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:40 SLIDINGWINDOW:4:10
		rm -r ${name}/${name}_1.fastq* -2 ${name}/${name}_2.fastq*
		pigz -p $SLURM_CPUS_PER_TASK ${name}/*.fq
		fi

	#cat ${name}/unpaired1.fq.gz ${name}/unpaired2.fq.gz > ${name}/unpaired.fq.gz
	#if [[ $(zcat unpaired.fq.gz | wc -l) > 0 ]]
	#then
	#	tophat_readsin="${name}/paired1.fq.gz ${name}/paired2.fq.gz,${name}/unpaired.fq.gz"
	#else
	#	tophat_readsin="${name}/paired1.fq.gz ${name}/paired2.fq.gz"
	#fi
	hisat2 -p $SLURM_CPUS_PER_TASK -x /data/Sherlock_Lung/JohnMce/CHM13/chm13v2.0 -1 ${name}/paired1.fq.gz -2 ${name}/paired2.fq.gz -U ${name}/unpaired1.fq.gz,${name}/unpaired2.fq.gz | samtools fastq -f 4 -F 512,1024 -0 /dev/null -s ${name}/${name}_loners.fq -1 ${name}/${name}_paired1.fq -2 ${name}/${name}_paired2.fq
	fi
	rm -r ${name}/paired1.fq* ${name}/unpaired1.fq* ${name}/paired2.fq* ${name}/unpaired2.fq*
	pigz -p $SLURM_CPUS_PER_TASK ${name}/*.fq
fi

for conf in "${confidences[@]}"
do
  kraken2 \
  --threads $SLURM_CPUS_PER_TASK -db ${db} \
  --paired --gzip-compressed --confidence ${conf} \
  --output ${name}/${name}-kraken-conf${conf}-pair.txt \
  --report ${name}/${name}-kraken-report-conf${conf}-pair.txt \
  ${name}/${name}_paired1.fq.gz ${name}/${name}_paired2.fq.gz

  kraken2 \
  --threads $SLURM_CPUS_PER_TASK -db ${db} \
  --gzip-compressed --confidence ${conf} \
  --output ${name}/${name}-kraken-conf${conf}-unpair.txt \
  ${name}/${name}_loners.fq.gz

  python ~/KrakenTools-master/make_kreport.py -i <(cat ${name}/${name}-kraken-conf${conf}-pair.txt ${name}/${name}-kraken-conf${conf}-unpair.txt) -t ${db}/ktaxonomy.tsv -o ${name}/${name}-kraken-report-conf${conf}.txt
	# run bracken shotgun parameters
	if [[ $Ssource == "RNA" || $Ssource == "DNA" ]]
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
