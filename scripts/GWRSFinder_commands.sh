#! /bin/bash
TALLYMER='gt'
GENOME='../blast_index/cbs138.fasta'
NAME='cbs138'
MIN='10'
BLAST='../20mers_minocc10/repeated_sequences_10.fasta'
ALLSTRAINS='../all_strains_blast_index/all_genomes.fasta'
outfmt1=(-outfmt '6 qseqid qstart qend sstart send nident mismatch gaps sstrand')
outfmt2=(-outfmt '6 qseqid sseqid length pident mismatch gaps')
### Index
	${TALLYMER} \
	suffixerator \
	-dna \
	-pl \
	-tis \
	-suf \
	-lcp \
	-v \
	-db ${GENOME} \
	-indexname ${NAME}

### Query
	${TALLYMER} \
	tallymer mkindex -scan -mersize 20 \
	-minocc ${MIN} \
	-esa ${NAME} \
	> ../20mers_minocc${MIN}/20mers_minocc_${MIN}.txt

### Index 
	${TALLYMER} \
	tallymer mkindex \
	-scan -mersize 20 \
	-minocc ${MIN} \
	-indexname tyr-${NAME} \
	-counts -pl \
	-esa ${NAME}
	
### 20mer position
	${TALLYMER} tallymer search -output qseqnum qpos counts sequence -tyr \
	tyr-${NAME} -q ${GENOME}  > ../20mers_minocc${MIN}/20mers_minocc${MIN}_pos.txt 

### add 20 to start position to get end position
	tallymer2bed.py \
	../20mers_minocc${MIN}/20mers_minocc${MIN}_pos.txt \
	../20mers_minocc${MIN}/20mers_minocc${MIN}_pos.bed \

### bedtools merge: merging only those that overlap or are bookended and are on the same strand
	bedtools merge \
	-i ../20mers_minocc${MIN}/20mers_minocc${MIN}_pos.bed -s -c 6 -o distinct \
	> ../20mers_minocc${MIN}/20mers_minocc${MIN}_merged.txt

### repeat_seq_3: make fasta file out of repeated sequences
	repeat_seq_3.py \
	$GENOME \
	../20mers_minocc${MIN}/20mers_minocc${MIN}_merged.txt \
	../20mers_minocc${MIN}/repeated_sequences_${MIN}.fasta

### blast against own database
	blastn -db ${GENOME} -query ${BLAST} "${outfmt1[@]}" > ../20mers_minocc${MIN}/unfiltered_blast.txt
	
### filter BLAST result
	repeatseq_blastfilter.py ../20mers_minocc${MIN}/unfiltered_blast.txt ../20mers_minocc${MIN}/filtered_blast.txt
	
### get sequence occurrence count
	count_blastrepeats.py ../20mers_minocc${MIN}/filtered_blast.txt ../20mers_minocc${MIN}/repeated_sequences_${MIN}.txt

### blast against database of all strains
	blastn -db ${ALLSTRAINS} -query ${BLAST} "${outfmt2[@]}" > ../20mers_minocc${MIN}/unfilteredblast_all_strains.txt
	
### filter blast result
	repeatseq_blastperfect.py ../20mers_minocc${MIN}/unfilteredblast_all_strains.txt \
	../20mers_minocc${MIN}/filteredblast_all_strains.txt \
	
### add filtered all strain blast result to repeated sequences file
	repeatseq_final.py ../20mers_minocc${MIN}/repeated_sequences_${MIN}.txt \
	../20mers_minocc${MIN}/filteredblast_all_strains.txt ../20mers_minocc${MIN}/repeated_sequencesfinal_${MIN}.txt


