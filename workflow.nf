#!/usr/bin/env nextflow

//params.references = "$projectDir/ncbi_tropheryma_whipplei/*.{fasta,fa,fna}"
params.references = "$projectDir/tests/*.{fasta,fa,fna}"
params.outdir = "$projectDir/results"
params.min_occurrence = 2

log.info """\
    GWRSFINDER   P I P E L I N E
    ===================================
    reference folder : ${params.references}
    outdir           : ${params.outdir}
    min occurrence       : ${params.min_occurrence}
    """
    .stripIndent()

process PROCESSGENOMES {
    label 'gwrsfinder'
    """ 
    simple python script that is run on each reference. 
    look for all fna files, get the name before '_', add to the fasta header. 
    """
    
    input:
    path fastafile

    output:
    path 'results.fasta'

    script:
    """
    python3 $projectDir/bin/Format_Fasta.py -i $fastafile -o results.fasta
    """
}

process CREATEDB {
    label 'ncbi'
    """
    Create a database containing all the genomes. 
    """
    input:
    path (genome_files, stageAs: "?/*") 

    output:
    tuple path("genomes.fasta"), path("genomes.fasta.*")

    script:
    """
    cat $genome_files > genomes.fasta
    makeblastdb -in genomes.fasta -dbtype nucl -parse_seqids
    """
}

process RUNTALLYMER {
    label 'gwrsfinder'
    """
    Create the tallymer indices and the kmer files.
    """

    input:
    path fastafile

    output:
    path "20mers_${params.min_occurrence}.txt"
    
    script:
    """
    gt suffixerator -dna -pl -tis -suf -lcp -v -db $fastafile -indexname $fastafile

    gt tallymer mkindex -scan -mersize 20 -minocc $params.min_occurrence -esa $fastafile

    gt tallymer mkindex -scan -mersize 20 -minocc $params.min_occurrence -indexname tyr-$fastafile -counts -pl -esa $fastafile

    gt tallymer search -output qseqnum qpos counts sequence -tyr tyr-$fastafile -q $fastafile  > 20mers_${params.min_occurrence}.txt
    """
}

process TALLYMER2BED {
    label 'gwrsfinder'
    """
    Convert the kmer file into a format suitable for bedtools.
    """

    input:
    path repeat_file

    output:
    path 'bed_formatted.txt'

    script:
    """
    python3 $projectDir/bin/Tallymer_2_Bed.py -i $repeat_file -o bed_formatted.txt -k ${params.min_occurrence}
    """
}

process BEDTOOLSMERGE {
    label 'gwrsfinder'
    """
    Merge overlapping or adjacent kmers with bedtools merge.
    """
    input:
    path bed_file
    
    output:
    path "20mers_minocc${params.min_occurrence}_merged.txt"

    script: 
    """
	bedtools merge -i $bed_file -c 6 -o distinct > 20mers_minocc${params.min_occurrence}_merged.txt
    """
}

process SEQJSON {
    label 'gwrsfinder'
    """
    Convert the merge file into a fasta file based on the sequences found in the original genome at 
    the given positions. 
    """

    input:
    path merged_bed_file
    path genome_file
    
    output:
    path "repeated_sequences_${params.min_occurrence}.json"
    
    script: 
    """ 
	python3 $projectDir/bin/MultiCopy_Parser.py -g $genome_file -i $merged_bed_file -o repeated_sequences_${params.min_occurrence}.json
    """
}

process CONSOLIDATESEQS {
    label 'gwrsfinder'
    publishDir params.outdir, mode:'copy'
    """
    consolidate all the repeated seqs from all the json files. 
    """
    input:
    // done because files are all the same name:https://stackoverflow.com/questions/73660749/nextflow-name-collision
    path(json_files, stageAs: "?/*") 

    output:
    path 'full_seqs.fasta'

    script:
    """
    python3 $projectDir/bin/JSON_2_Fasta.py $json_files
    """
}

process BLASTSEQS {
    label 'ncbi'
    
    """
    Blast the multi-copy fasta file against the genomic database. 
    """
    input:
    path multicopy_fasta_file
    tuple path(db_file), path(index_files)

    output:
    path 'unfiltered_blast.txt'

    script:
    """
    blastn -db ${db_file[0]} -query $multicopy_fasta_file -outfmt '6 qseqid sseqid length pident mismatch gaps' \
    > unfiltered_blast.txt
    """

}

process AGGREGATERESULTS {
    label 'gwrsfinder'
    publishDir params.outdir, mode:'copy'
    """
    From the results of the blast, count the occurrences of each multicopy sequence in each genome and 
    then print out to a csv. 
    """
    input:
    path multicopy_fasta_file
    path unfiltered_results

    output:
    path 'pcr_targets.json'

    script:
    """
    python3 $projectDir/bin/Aggregate_Blast.py -m $multicopy_fasta_file -b $unfiltered_results -o 'pcr_targets.json'
    """
}

workflow {

    Channel
        .fromPath(params.references, checkIfExists: true)
        .set { reference_files }
  
    reference_ch = PROCESSGENOMES(reference_files)

    db_file = CREATEDB(reference_ch.collect())
    repeat_file = RUNTALLYMER(reference_ch)
    bed_compatible_repeat_file = TALLYMER2BED(repeat_file)
    merged_file = BEDTOOLSMERGE(bed_compatible_repeat_file)
    
    json_files = SEQJSON(merged_file, reference_ch)
    multicopy_fasta_file = CONSOLIDATESEQS(json_files.collect())
    unfiltered_results = BLASTSEQS(multicopy_fasta_file, db_file)
    AGGREGATERESULTS(multicopy_fasta_file, unfiltered_results)
}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone!" : "Oops .. something went wrong" )
}