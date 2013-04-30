#!/usr/bin/env nextflow

/*
 * Defines the pipeline parameters.
 * The values in the 'params' map can be overridden on the command line by specifying a
 * option prefixed with a double '-' char, for example
 *
 * $ nextflow piper.nf --query=<path to your file name>
 *
 */

params['query-chunk-len'] = 100
params['query'] = "${HOME}/workspace/piper/tutorial/5_RNA_queries.fa"
params['genomes-db'] = "${HOME}/workspace/piper/tutorial/db"
params['max-threads'] = Runtime.getRuntime().availableProcessors()

// these parameters are mutually exclusive
// Input genome can be specified by
// - genomes-file: a file containing the list of genomes FASTA to be processed
// - genomes-list: a comma separated list of genomes FASTA file
// - genomes-folder: a directory containing a folder for each genome FASTA file
params['genomes-file'] = null
params['genomes-list'] = null
params['genomes-folder'] = "${HOME}/workspace/piper/tutorial/genomes/"


queryFile = file(params.query)
dbpath = file(params['genomes-db'])

if( !dbpath.exists() ) {
    log.warn "Creating genomes-db path: $dbpath"
    if( !dbpath.mkdirs() ) {
        exit 1, "Cannot create genomes-db path: $dbpath -- check file system permissions"
    }
}


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


allGenomes = [:]

// -- when the provided source path is a FILE
//      each line represent the path to a genome file
if( params['genomes-file'] ) {
    def sourcePath = new File(params['genomes-file'])
    if( sourcePath.isEmpty() ) {
        exit 1, "Not a valid input genomes descriptor file: ${sourcePath}"
    }

    // parse the genomes input file files (genome-id, path to genome file)
    int count=0
    sourcePath.eachLine { line ->
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
            exit 2, "Missing input genome file: $fasta"
        }

        allGenomes[ genomeId ] = [
                genome_fa: fasta,
                chr_db: new File(dbpath,"${genomeId}/chr"),
                blast_db: new File(dbpath, "${genomeId}/blastdb")
            ]
    }
}
else if( params['genomes-list'] ) {

    def count=0
    def files = params['genomes-list'].split(',').collect { new File(it.trim()) }

    files.each { fasta ->
        if( !fasta.exists() ) {
            exit 3, "Missing genome file: $fasta"
        }

        def genomeId = "gen${++count}"
        allGenomes[ genomeId ] = [
                genome_fa: fasta,
                chr_db: new File(dbpath,"${genomeId}/chr"),
                blast_db: new File(dbpath, "${genomeId}/blastdb")
            ]

    }

}
else if( params['genomes-folder'] ) {
    def sourcePath = new File(params['genomes-folder'])
    if( !sourcePath.exists() || sourcePath.isEmpty() ) {
        exit 4, "Not a valid input genomes folder: ${sourcePath}"
    }

    sourcePath.eachDir { File path ->
        def fasta = path.listFiles().find{ File file -> file.name.endsWith('.fa') }
        if( fasta ) {
            log.info "Processing => ${path.name} - ${fasta}"
            allGenomes[ path.name ] = [
                    genome_fa: fasta,
                    chr_db: new File(dbpath,"${path.name}/chr"),
                    blast_db: new File(dbpath, "${path.name}/blastdb")
                ]
        }
    }
}

else {
    exit 5, "No input genome(s) provided -- Use one of the following CLI options 'genomes-file' or 'genomes-list' or 'genomes-folder' "
}

if( !allGenomes ) {
    exit 6, "No genomes found in path"
}

// get all genomes ID found and put into a list
formatName = allGenomes.keySet()


/*
 * Split the query input file in many small files (chunks).
 *
 * The number of sequences in each chunk is controlled by the parameter 'query-chunk-len'
 * The chunk files are saved in a local folder define by the variable 'querySplits'
 *
 */

// create a folder that may be cached, using the 'queryFile' and the number chunks as cache key
querySplits = cacheableDir([queryFile, params.'query-chunk-len'])

if( querySplits.isEmpty() ) {
    log.info "Splitting query file: $queryFile .."
    chunkCount=0
    queryFile.chunkFasta( params.'query-chunk-len' ) { sequence ->
        def file = new File(querySplits, "seq_${chunkCount++}")
        file.text = sequence
    }
    log.info "Created $chunkCount input chunks to path: ${querySplits}"
}
else {
    log.info "Cached query splits > ${querySplits.list().size()} input query chunks"
}



allQueryIDs = []
queryFile.chunkFasta() { String chunk ->
    allQueryIDs << chunk.readLines()[0].substring(1)
}


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
    threads params['max-threads']

    """
    set -e
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
        SEQID=`grep -o -E "^>\\S+" $x | tr -d ">" | sed 's/[\\>\\<\\/\\''\\:\\\\]/_/'`
        mv $x ${CHR_DB}/$SEQID;
        done

    fi

    echo $NAME > blastName

    """
}



/*
 * Iterate over the query chunks and create a pair (genome name, chunk file) for each of them
 */
blastId = channel()
blastQuery = channel()

blastName.each {

    def name = it.text.trim()
    querySplits.eachFile { chunk ->
        log.info "Blasting > $name - chunk: $chunk"
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
    input blastId
    input blastQuery
    output exonerateId
    output exonerateQuery
    output blastResult
    threads params['max-threads']

    """
    set -e
    echo ${blastId} > exonerateId
    blastn -db ${allGenomes[blastId].blast_db}/db -query ${blastQuery} -outfmt '$outFmt' > blastResult
    ln -s ${blastQuery} exonerateQuery
    """

}

/*
 * Given the blast output execute the 'exonerate' step
 */

exonerateOut = channel()

task ('exonerate') {
    input exonerateId
    input exonerateQuery
    input blastResult
    output '*.fa': exonerateOut
    threads params['max-threads']

    """
    specie='${exonerateId.text.trim()}'
    chr=${allGenomes[exonerateId.text.trim()].chr_db}
    ## apply exonerate
    exonerateRemapping.pl -query ${exonerateQuery} -mf2 $blastResult -targetGenomeFolder $chr -exonerate_lines_mode 1000 -exonerate_success_mode 1 -ner no

    ## exonerateRemapping create a file named 'blastResult.fa'
    ## split the exonerate result into single files
    ${split_cmd} blastResult.fa '%^>%' '/^>/' '{*}' -f .seq_ -n 5
    mv blastResult.fa .blastResult.fa

    ## rename the seq_xxx files so that the file name match the seq fasta id
    ## plus append the specie to th sequence id
    for x in .seq_*; do
      SEQID=`grep '>' \$x`
      FILENAME=`grep '>' \$x | sed 's/^>\\(.*\\)_hit\\d*.*\$/\\1/'`
      printf "\${SEQID}_${specie}\\n" > \${FILENAME}.fa
      cat \$x | grep -v '>' >> \${FILENAME}.fa
    done
    """
}


fastaToMerge = new Channel()
exonerateOut.filter { file -> file.baseName in allQueryIDs  } .into (fastaToMerge)

fastaToAlign = merge('prepare_mfa') {

    input fastaToMerge
    output '*.mfa'

    """
    # Extract the file name w/o the extension
    fileName=\$(basename "$fastaToMerge")
    baseName="\${fileName%.*}"

    # Only the first time append the query sequence
    if [ ! -e \$baseName.mfa ]; then
    perl -n -e '$on=(/^>('\$baseName')\$/) if (/^>/); print $_ if ($on);' $queryFile > \$baseName.mfa
    fi

    # Append the exonerate result
    cat $fastaToMerge >> \$baseName.mfa
    """
}

alignment = task('align') {
    input fastaToAlign
    output '*.aln'

    """
    t_coffee -in $fastaToAlign -method slow_pair -n_core 1
    """
}

similarity = merge('similarity') {
    input alignment
    output '*'

    """
    fileName=\$(basename "$alignment")
    baseName="\${fileName%.*}"
    t_coffee -other_pg seq_reformat -in $alignment -output sim > \$baseName
    """
}

simFolder = val()
similarity.whenBound { file -> simFolder << file.parent }

task ('matrix') {
    echo true
    input simFolder

    """
    sim2matrix.pl -query $queryFile -data_dir $simFolder -genomes_dir $dbpath
    """
}
