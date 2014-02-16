genomes = Channel.path( params.genomes )

  process formatBlastDatabases {

    storeDir 'blastdb'

    input:
    file specie from genomes

    output:
    file "${dbName}.*" into blastDb

    script:
    dbName = specie.baseName
    """
    makeblastdb -dbtype prot -in ${specie} -out ${dbName}
    """

  }

blastDb.subscribe { println it  }

