# GWRSFinder PIPELINE
    Purpose: To interrogate small genomes and identify conserved multicopy sequences for PCR targeting. For a given species, a reference genome is selected and the software Tallymer (part of the genometools software, v1.5.9) identifies kmers (20mers are used in this pipeline) repeated within the reference genome and their position in the genome. A custom python script (tallymer2bed.py) generates a bed file (20mers_minoccX_pos.bed) from the tallymer output and Bedtools merge (v2.26.0) merges overlapping or bookended 20mers and identifies positions of merged sequences in the reference genome, the output of which is saved to a text file (20mers_minoccX_merged.txt). From the positions of the multicopy sequences given by the bedtools file, a custom Python script is used to retrieve the sequences from the genome and create a fasta file containing these sequences (repeat_seq_3.py, repeated_sequences_X.fasta). BLAST command line tools (v2.6.0) is used to create a local database for the reference genome and available genomes of other strains from the same species and is used to search for the multicopy sequences in this database. Custom Python scripts (repeatseq_blastfilter.py, count_blastrepeats.py, repeatseq_blastperfect.py, repeatseq_final.py) use the output from the BLAST search to generate a final list of multicopy sequences, their copy number in the reference genome and their copy number in the genomes of other strains from the same species (repeated_sequencesfinal_X.txt)

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
4. Add scripts in ./script to PATH and ensure scripts are executable

--------------------------------------------------------------------------------------------------
# INSTRUCTIONS
    Examples Candida glabrata and Mycobacterium tuberculosis are in the in test directory. 
    The steps below contain command line instructions that must be run and the directory they must 
    be run in (dir:). 
--------------------------------------------------------------------------------------------------
# folder structure for pipeline:
    see ./test/
        ./blast_index #contains the reference genome
        ./all_strains_blast_index #contains all other genomes to be checked for the multi-copy sequences
        ./tallymer_index #will contain tallymer files
        ./20mers_minocc$X #will contain output files

# If genomes are multiple fasta sequences, need to change to single fasta sequence using the following commands:
1. Convert multifasta files to single fasta files for all genomes (reference and non-reference)
    
    C. glabrata 
    
    dir: ./blast_index:

        $ single_fasta.py cbs138.fna cbs138

    dir: ./all_strains_blast_index

        $ single_fasta.py w10d4.fna w10d4
        $ single_fasta.py dsy565.fna dsy565
        $ single_fasta.py dsy562.fna dsy562

    This is not necessary for M. tuberculosis as genomes are a single fasta

2. Put all (non-reference) genomes into one file
   
    dir: ./all_strains_blast_index
   
    C. glabrata:
        
         $ cat w10d4.fasta dsy565.fasta dsy562.fasta >> all_genomes.fasta
    M. tuberculosis:

   	    $ cat GCF_003287165.1.fna GCF_001544705.1.fna GCF_003287145.1.fna >> all_genomes.fasta

# Create separate BLAST databases for the reference genome and the other genomes:
1. dir: ./blast_index
    
    C. glabrata:

        $ makeblastdb -in cbs138.fasta -dbtype nucl -parse_seqids
    M. tuberculosis:
        
        $ makeblastdb -in GCF_000195955.2.fna -dbtype nucl -parse_seqids
2. dir: ./all_strains_blast_index
    
    C. glabrata and M. tuberculosis:
	    
        $ makeblastdb -in all_genomes.fasta -dbtype nucl -parse_seqids

# Remaining command line instructions are run from the shell script 'GWRSFinder_commands.sh'
1. dir: ./tallymer_index
    
    change the following variables in GWRSFinder_commands.sh to fit the genome(s) being searched: 
        
        GENOME='../cbs138.fna' #reference genome file
        NAME='cbs138' #reference genome name
        MIN='10' #minimum threshold the 20mer must be present in the genome to be picked up
        BLAST='../20mers_minocc10/repeated_sequences_10.fasta' #file where blast result is saved
    
    example: C. glabrata:	
	    
        $ TALLYMER='gt'
	    $ GENOME='../blast_index/cbs138.fasta'
	    $ NAME='cbs138'
	    $ MIN='10'
	    $ BLAST='../20mers_minocc10/repeated_sequences_10.fasta'
	    $ ALLSTRAINS='../all_strains_blast_index/all_genomes.fasta'
	    $ outfmt1=(-outfmt '6 qseqid qstart qend sstart send nident mismatch gaps sstrand')
	    $ outfmt2=(-outfmt '6 qseqid sseqid length pident mismatch gaps')

    
    example: M. tuberculosis:
	
        $ TALLYMER='gt'
	    $ GENOME='../blast_index/GCF_000195955.2.fna'
	    $ NAME='GCF_000195955.2'
	    $ MIN='5'
	    $ BLAST='../20mers_minocc5/repeated_sequences_5.fasta'
	    $ ALLSTRAINS='../all_strains_blast_index/all_genomes.fasta'
	    $ outfmt1=(-outfmt '6 qseqid qstart qend sstart send nident mismatch gaps sstrand')
	    $ outfmt2=(-outfmt '6 qseqid sseqid length pident mismatch gaps')
	

--------------------------------------------------------------------------------------------------
# RESULT
--------------------------------------------------------------------------------------------------
The result will be in the the 'repeated_sequencesfinal_X.txt' file, where X refers to the minimum
number of repeats chosen. This is a tab delimited file containing:

'sequence'    'sequence length' 'copy number' 'copy number in non-reference genomes'

'copy number' refers to copy number in reference genome.

--------------------------------------------------------------------------------------------------
# Limitations to consider
--------------------------------------------------------------------------------------------------
Depending on the abundance (or lack thereof) of multicopy sequences within a given reference genome,
the minimum copy number parameter of the 20mer should be adjusted. For instance, if the highest
number of times a 20mer is present in a genome is 5, this will not show up if the minimum number is 
set to 10. Additionally, this minimum number relates to the minimum times the 20mer must be present, 
not the final multicopy sequences. For this reason, the number of times a final multi-copy sequence 
is be present may be much lower than a given 20mer within that final sequence.

Secondly, the final tab delimited file containing the mult-copy sequence list may have many sequences it,
and many of these sequences may be very similar give or take a few nucleotides. For example, a sequence 
that is 60 nucleotides in length may be present 10 times, while another sequence may contain the exact 
same 60 nucleotides plus an additional 20 nucleotides (80 nucleotides in length total) but only be present
in the genome 5 times. Both of these sequences may be in the final list and manual parsing through the list
is required to check this. It is kept this way for primer design. The 60 nucleotide sequence is more desirable
because it is of higher copy number, however it may not be possible to design primers on then and so the 80 
length sequence might be the better option. 

Lastly, this pipeline has only been tested on bacterial and fungal genomes. It is not recommended to try anything 
larger