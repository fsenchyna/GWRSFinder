#!/usr/bin/env nextflow

log.info """\
    GWRSFINDER   P I P E L I N E
    ===================================
    reference folder : ${params.references}
    outdir           : ${params.outdir}
    min occurrence       : ${params.min_occurrence}
    accession format    : ${params.accession_id}
    """
    .stripIndent()

process PARSEACCESSION {
    """ Parse the accession id from the fasta file name. """

    tag "Extracting accesssion id from $fastafile "
    label 'gwrsfinder'
    
    input:
    path fastafile   

    output:
    stdout

    script:
    """
    re="(${params.accession_id}[^_]*?)[_.].*[fasta|fa|fna]"
    if [[ $fastafile =~ \$re ]]; then echo -n \${BASH_REMATCH[1]}; fi
    """
}

process PROCESSGENOMES {
    """ Concatenate accession id to fasta headers. """

    tag "Concatenating $accession_id to $fastafile headers. "
    label 'gwrsfinder'

    
    input:
    path fastafile 
    val accession_id

    output:
    path "${accession_id}.fasta"

    script:
    """
    python3 $projectDir/bin/FastaHeaderFormatter/FastaHeaderFormatter.py -i $fastafile -id $accession_id \
    -o ${accession_id}.fasta
    """
}

process CREATEDB {
    """ Create a database for blast containing all the reference files. """

    tag "Creating the blast database. "
    label 'ncbi'

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
    """ Create the tallymer indices and the kmer files. """

    tag "Creating the tallymer index and kmer files "
    label 'gwrsfinder'

    input:
    path fastafile

    output:
    path "20k_${params.min_occurrence}m.txt"
    
    script:
    """
    gt suffixerator -dna -pl -tis -suf -lcp -v -db $fastafile -indexname $fastafile

    gt tallymer mkindex -scan -mersize 20 -minocc $params.min_occurrence -esa $fastafile

    gt tallymer mkindex -scan -mersize 20 -minocc $params.min_occurrence -indexname tyr-$fastafile -counts -pl -esa $fastafile

    gt tallymer search -output qseqnum qpos counts sequence -tyr tyr-$fastafile -q $fastafile  > 20k_${params.min_occurrence}m.txt
    """
}

process TALLYMER2BED {
    """ Convert kmer file into a format suitable for bedtools merge. """

    tag "converting $repeat_file kmer file to bed format. "   
    label 'gwrsfinder'

    input:
    path repeat_file

    output:
    path "bed_${repeat_file}"

    script:
    """
    python3 $projectDir/bin/TallymerToBed/TallymerToBed.py -i $repeat_file -o bed_${repeat_file} -k ${params.min_occurrence}
    """
}

process BEDTOOLSMERGE {
    """ Merge overlapping or adjacent kmers with bedtools merge. """

    tag " Running bedtools merge on $bed_file "  
    label 'gwrsfinder'

    input:
    path bed_file
    
    output:
    path "merged_$bed_file"

    script: 
    """
	bedtools merge -i $bed_file -c 6 -o distinct > merged_$bed_file
    """
}

process SEQJSON {
    """ Extract sequences from genome based on bed file. """

    tag " Extracting multicopy sequences from $genome_file"  
    label 'gwrsfinder'

    input:
    path merged_bed_file
    path genome_file
    val accession_id
    
    output:
    path "${accession_id}.json"
    
    script: 
    """ 
	python3 $projectDir/bin/MultiCopyParser/MultiCopyParser.py -g $genome_file -i $merged_bed_file \
    -o ${accession_id}.json
    """
}

process CONSOLIDATESEQS {
    """ Consolidate multicopy sequences across all genomes. """

    tag " Consolidating multicopy sequences across all genomes. "  
    label 'gwrsfinder'
    publishDir params.outdir, mode:'copy'

    input:
    // done because files are all the same name:https://stackoverflow.com/questions/73660749/nextflow-name-collision
    path json_files

    output:
    path 'full_seqs.fasta'

    script:
    """
    python3 $projectDir/bin/JsonToFasta/JsonToFasta.py -i $json_files -o full_seqs.fasta
    """
}

process BLASTSEQS {
    """ Blast multi-copy sequences against the genomic database. """

    tag " Blasting multi-copy sequences against the genomic database. "  
    label 'ncbi'

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
    """ Consolidate the results into a json file. """

    tag " Consolidating the results. " 
    label 'gwrsfinder'
    publishDir params.outdir, mode:'copy'

    input:
    path multicopy_fasta_file
    path unfiltered_results

    output:
    path 'pcr_targets.json'

    script:
    """
    python3 $projectDir/bin/TransformBlast/TransformBlast.py -b $unfiltered_results -o pcr_targets.json
    """
}

workflow {

    Channel
        .fromPath(params.references, checkIfExists: true)
        .set { reference_files }

    accession_ids = PARSEACCESSION(reference_files)
    accession_ids.view { "Accession: ${it}" }

    reference_ch = PROCESSGENOMES(reference_files, accession_ids)

    db_file = CREATEDB(reference_ch.collect())
    repeat_file = RUNTALLYMER(reference_ch)
    bed_compatible_repeat_file = TALLYMER2BED(repeat_file)
    merged_file = BEDTOOLSMERGE(bed_compatible_repeat_file)
    
    json_files = SEQJSON(merged_file, reference_ch, accession_ids)
    multicopy_fasta_file = CONSOLIDATESEQS(json_files.collect())
    unfiltered_results = BLASTSEQS(multicopy_fasta_file, db_file)
    AGGREGATERESULTS(multicopy_fasta_file, unfiltered_results)
}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone!" : "Oops .. something went wrong" )
}