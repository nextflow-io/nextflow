genomes = Channel.path( params.genomes )

  process formatBlastDatabases {

    storeDir 'blastdb'

    input:
    file species from genomes

    output:
    file "${dbName}.*" into blastDb

    script:
    dbName = species.baseName
    """
    makeblastdb -dbtype prot -in ${species} -out ${dbName}
    """

  }

blastDb.subscribe { println it  }

