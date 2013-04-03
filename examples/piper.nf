#!/bin/env nextflow

echo true

/*
 * Defines the pipeline parameters.
 * The values in the 'params' map can be overridden on the command line by specifying a
 * option prefixed with a double '-' char, for example
 *
 * $ nextflow piper.nf --query=<path to your file name>
 */
params['query-chunk-len'] = 100
params['query'] = "${HOME}/workspace/piper/tutorial/5_RNA_queries.fa"
params['genomes'] = "${HOME}/workspace/piper/tutorial/genomes/"
params['genomes-db'] = "${HOME}/workspace/piper/tutorial"

queryFile = new File( params.query )
querySplits = new File('./splits')
dbpath = new File(params.'genomes-db')
sourceGenomesPath = new File(params.genomes)


/*
 * Find out all the genomes files in the specified directory.
 *
 * More in detail teh 'sourceGenomesPath' points to a directory having a
 * sub-folder for each genome it is required to process.
 *
 * Each sub-folder must contain the genome FASTA file to be processed.
 *
 * The sub-folder name is used to identify the genome in the computation.
 *
 * All the genomes names found in this path are put in a list named 'formatName',
 * which control the pipeline execution.
 *
 */


formatName = []
allGenomes = [:]

// -- when the provided source path is a FILE
//      each line represent the path to a genome file
if( sourceGenomesPath.isFile() ) {
    // parse the genomes input file files (genome-id, path to genome file)
    int count=0
    sourceGenomesPath.eachLine { line ->
        count++
        def genomeId
        def path

        def items = line.trim().split(/\s/)
        if( items.size() > 1 ) {
            (genomeId, path) = items
        }
        else if( items.size() ==1  ){
            genomeId = "gen${count}"
            path = items[0]
        }
        else {
            return
        }

        def fasta = new File(path)
        if( !fasta.exists() ) {
            println "Missing input genome file: $fasta"
            exit 1
        }

        allGenomes[ genomeId ] = [
                genome_fa: fasta,
                chr_db: new File(dbpath,"${genomeId}/chr"),
                blast_db: new File(dbpath, "${genomeId}/blastdb")
            ]
        formatName << genomeId
    }

}
else if( sourceGenomesPath.isDirectory() ) {

    sourceGenomesPath.eachDir { File path ->
        def fasta = path.listFiles().find{ File file -> file.name.endsWith('.fa') }
        if( fasta ) {
            println "Processing => ${path.name} - ${fasta}"
            allGenomes[ path.name ] = [
                    genome_fa: fasta,
                    chr_db: new File(dbpath,"${path.name}/chr"),
                    blast_db: new File(dbpath, "${path.name}/blastdb")
                ]
            formatName << path.name
        }
    }
}
else {
    println "Invalid genomes source file: ${sourceGenomesPath}"
    exit 1
}


if( !allGenomes ) {
    println "No genomes found in path: ${sourceGenomesPath}"
    exit 1
}


/*
 * Split the query input file in many small files (chunks).
 *
 * The number of sequences in each chunk is controlled by the parameter 'query-chunk-len'
 * The chunk files are saved in a local folder define by the variable 'querySplits'
 *
 */
println "Splitting query file: $queryFile .."
querySplits.with {
    if( exists() && !isEmpty() ) { deleteDir() }
    if( !exists() && !mkdirs() ) {
        println "Unable to create query splits folder: $querySplits"
        exit 2
    }
}

int chunkCount=0
queryFile.chunkFasta( params.'query-chunk-len' ) { sequence ->
    def file = new File(querySplits, "seq_${chunkCount++}")
    file.text = sequence
}
println "Created $chunkCount input chunks to path: ${querySplits}"


/*
 * Create the required databases (BLAST,CHR) if they does not exists.
 *
 * This task is executed for each genome in the list 'formatName'
 * The tasks 'sends' out the name of the genome to be processed
 * by the next step in the pipeline using the variable 'blastName'
 */


def split_cmd = (System.properties['os.name'] == 'Mac OS X' ? 'gcsplit' : 'csplit')

task('format') {
    input formatName
    output blastName
    threads 8

    """
    NAME=${formatName}
    FASTA=${allGenomes[formatName].genome_fa}
    CHR_DB=${allGenomes[formatName].chr_db}
    BLAST_DB=${allGenomes[formatName].blast_db}

    ## Create the BLAST db if they does not exist
    if [[ ! `ls -A ${BLAST_DB} 2>/dev/null` ]]; then

        ## Create the target folder
        mkdir -p ${BLAST_DB}

        ## Format the BLAST DB
        makeblastdb -dbtype nucl -in ${FASTA} -out ${BLAST_DB}/db
    fi


    ## Create the CHR database if does not exist
    if [[ ! `ls -A ${CHR_DB} 2>/dev/null` ]]; then

        ## split the fasta in a file for each sequence 'seq_*'
        ${split_cmd} ${FASTA} '%^>%' '/^>/' '{*}' -f seq_ -n 5

        ## create the target folder
        mkdir -p ${CHR_DB}

        ## rename and move to the target folder
        for x in seq_*; do
        SEQID=`grep -o -E "^>\\w+" $x | tr -d ">"`
        mv $x ${CHR_DB}/$SEQID;
        done

    fi

    echo $NAME > blastName

    """
}



/*
 * Iterate over the query chunks and create a pair (genome name, chunk file) for each of them
 */
blastId = queue()
blastQuery = queue()

blastName.each {

    def name = it.text.trim()
    querySplits.eachFile { chunk ->
        println "Blasting '$name' - chunk: $chunk"
        synchronized(this) {
            blastId << name
            blastQuery << chunk.absoluteFile
        }
    }

}



/*
 * Implements the BLAST step
 */

outFmt = '6 qseqid sseqid evalue score qgi bitscore length nident positive mismatch pident ppos qacc gaps gaopen qaccver qlen qframe qstart qend sframe sstart send'


task ('blast') {
    echo true
    input blastId
    input blastQuery
    output exonerateId
    output exonerateQuery
    output blastResult
    threads 8

    """
    echo ${blastId} > exonerateId
    blastn -db ${allGenomes[blastId].blast_db}/db -query ${blastQuery} -outfmt '$outFmt' > blastResult
    ln -s ${blastQuery} exonerateQuery
    """

}

/*
 * Given the blast output execute the 'exonerate' step
 */
task ('exonerate') {
    input exonerateId
    input exonerateQuery
    input blastResult
    threads 8

    """
    chr=${allGenomes[exonerateId.text.trim()].chr_db}
    exonerateRemapping.pl -query ${exonerateQuery} -mf2 $blastResult -targetGenomeFolder $chr -exonerate_lines_mode 1000 -exonerate_success_mode 1 -ner no
    """

}

