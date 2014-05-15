
proteins = Channel.fromPath( "examples/data/p?.fa" ).buffer(size:2, remainder: true)

process blastThemAll {
    echo true

    input:
    file x from proteins

    "echo $x"

}
