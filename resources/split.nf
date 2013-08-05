import nextflow.util.DxFile

file = new DxFile(name:'sample.fa', id:'file-B7vQg7j0j583xBb7QJV00284')
//file = file("$HOME/sample.fa")

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
    rev $chunks
    """
}

