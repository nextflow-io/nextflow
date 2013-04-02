#!/bin/env nextflow

DB = '/Users/ptommaso/tools/blast-db/pdb/pdb'
fileName = '/Users/ptommaso/sample.fa'

seq = queue()
new File(fileName).chunkFasta { seq << it }


task {

    input '-': seq
    output blastResult

    """
    cat - | blastp -db $DB -query - -outfmt 6 > blastResult
    """

}


merge {

    input blastResult
    output allBlast

        """
        cat ${blastResult} >> allBlast
        """

}

task {
    input allBlast

        """
        sort $allBlast
        """
}


println read(result)
