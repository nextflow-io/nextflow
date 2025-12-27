nextflow.preview.types = true

process foo {

    output:
    out1: Path = file("out1.txt")
    out2: Path = file("out2.txt")

    topic:
    file("version.txt") >> 'versions'

    script:
    """
    echo "version" > version.txt
    echo "result1" > out1.txt
    echo "result2" > out2.txt
    """
}

workflow {
    (out1, out2) = foo()

    channel.topic('versions')
    | mix( out1 )
    | mix( out2 )
    | map { f -> f.text.trim() }
    | collectFile(name: 'result.txt', sort: true, newLine: true, storeDir: '.')
}
