
proteins = Channel.path( "examples/data/p?.fa" ).buffer(size:2)

process blastThemAll {
    echo true

    input:
    file x from proteins

    "echo $x"

}
