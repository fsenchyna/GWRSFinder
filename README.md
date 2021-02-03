# GWRSFinder PIPELINE
    Find multi-copy sequences in a reference genome and search for those multi-copy sequences in a list of other genomes,
    typically those genomes are from different isolates of the same organism as the reference genome. 

--------------------------------------------------------------------------------------------------
# INSTALLATION REQUIREMENTS 
--------------------------------------------------------------------------------------------------
1. genometools (v.1.5.9) 
    Installation source: http://genometools.org/pub/
2. BLAST+ command line applicatons (v2.6.0) 
    Installation source: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/
    Installation instructions: https://www.ncbi.nlm.nih.gov/books/NBK279671/#introduction.MacOSX
3. Bedtools (v2.26.0)
    Installation source: https://github.com/arq5x/bedtools2/releases
4. Add scripts in ./script to PATH and make the scripts executables

--------------------------------------------------------------------------------------------------
# INSTRUCTIONS
    Instructions for examples Candida glabrata and Mycobacterium tuberculosis in test folder. Contains
    command line instructions that must be run and the folder they must be run in (dir:). 
--------------------------------------------------------------------------------------------------
# folder structure for pipeline:
see ./test
    ./blast_index contains the reference genome
    ./all_strains_blast_index contains all other genomes to be checked for the multi-copy sequences
    ./tallymer_index will contain tallymer files
    ./20mers_minocc$X will contain output files

# If genomes are in multi-fasta format, need to change to single fasta sequence using the following commands:
convert multifasta files to single fasta files for all genomes (reference and non-reference)
C. glabrata (M. tuberculosis not necessary, all files are single fasta), dir: ./blast_index:
    $  single_fasta.py cbs138.fna cbs138
dir: ./all_strains_blast_index
    $  single_fasta.py w10d4.fna w10d4
    $  single_fasta.py dsy565.fna dsy565
    $  single_fasta.py dsy562.fna dsy562
put all (non-reference) genomes into one file
dir: ./all_strains_blast_index
C. glabrata:
    $ cat w10d4.fasta dsy565.fasta dsy562.fasta >> all_genomes.fasta
M. tuberculosis
   	$ cat GCF_003287165.1.fna GCF_001544705.1.fna GCF_003287145.1.fna >> all_genomes.fasta

# Create separate BLAST databases for the reference genome and the other genomes:
dir: ./blast_index
C. glabrata:
    $ makeblastdb -in cbs138.fasta -dbtype nucl -parse_seqids
M. tuberculosis:
    $ makeblastdb -in GCF_000195955.2.fna -dbtype nucl -parse_seqids
dir: ./all_strains_blast_index
C. glabrata and M. tuberculosis:
	$ makeblastdb -in all_genomes.fasta -dbtype nucl -parse_seqids

# Remaining command line instructions are run from the shell script 'GWRSFinder_commands.sh'
dir: ./tallymer_index
change the following variables in GWRSFinder_commands.sh to fit the genome(s) being searched: 
    GENOME='../cbs138.fna' # reference genome file
    NAME='cbs138' # reference genome name
    MIN='10' # minimum threshold the 20mer must be present in the genome to be detected
    BLAST='../20mers_minocc10/repeated_sequences_10.fasta' # file where blast result is saved
    
example: C. glabrata:	
	TALLYMER='gt'
	GENOME='../blast_index/cbs138.fasta'
	NAME='cbs138'
	MIN='10'
	BLAST='../20mers_minocc10/repeated_sequences_10.fasta'
	ALLSTRAINS='../all_strains_blast_index/all_genomes.fasta'
	outfmt1=(-outfmt '6 qseqid qstart qend sstart send nident mismatch gaps sstrand')
	outfmt2=(-outfmt '6 qseqid sseqid length pident mismatch gaps')

    
example: M. tuberculosis:
	TALLYMER='gt'
	GENOME='../blast_index/GCF_000195955.2.fna'
	NAME='GCF_000195955.2'
	MIN='5'
	BLAST='../20mers_minocc5/repeated_sequences_5.fasta'
	ALLSTRAINS='../all_strains_blast_index/all_genomes.fasta'
	outfmt1=(-outfmt '6 qseqid qstart qend sstart send nident mismatch gaps sstrand')
	outfmt2=(-outfmt '6 qseqid sseqid length pident mismatch gaps')
	

--------------------------------------------------------------------------------------------------
# RESULT
--------------------------------------------------------------------------------------------------
The result will be in the the 'repeated_sequencesfinal_X.txt' file, where X refers to the minimum
number of repeats chosen. This is a tab delimited file containing:

'sequence'    'sequence length' 'copy number' 'copy number in non-reference genomes'

'copy number' refers to copy number in reference genome.