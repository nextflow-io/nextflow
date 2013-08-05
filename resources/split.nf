import nextflow.util.DxFile

file = new DxFile(name:'sample.fa', id:'file-B7zq7j80FqXB4J080VV0023Y')
//file = file("$HOME/vagrant_machine_precise64/nextflow/sample.fa")

chunks = channel()

task {

    input file
    output 'seq_*': chunks

    """
    csplit $file '%^>%' '/^>/' '{*}' -f seq_
    """

}


task {
    input chunks
    echo true

    """
    echo "CHUNKS >> ${chunks}"
    rev $chunks
    """
}

