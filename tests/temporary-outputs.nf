
process foo {
    input:
    val meta_id
    output:
    tuple val(meta_id), path('a.txt', temporary: true)

    script:
    """
    touch a.txt
    echo 'foo was here' >> a.txt
    """
}

process bar {
    input:
    tuple val(meta_id), path('a.txt')
    output:
    tuple val(meta_id), path('b.txt', temporary: true)

    script:
    """
    cat a.txt > b.txt
    echo 'bar was here' >> b.txt
    """
}

process baz {
    publishDir '.'

    input:
    tuple val(meta_id), path('a.txt'), path('b.txt')
    output:
    tuple val(meta_id), path('c.txt')

    script:
    """
    cat b.txt > c.txt
    echo 'baz was here' >> c.txt
    """
}

workflow {
    meta_ids = Channel.of( '1', '2', '3' )
    ch_a = foo(meta_ids)
    ch_b = bar(ch_a)
    ch_ab = ch_a.join(ch_b)
    ch_c = baz(ch_ab)
}
