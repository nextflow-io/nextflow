params.str = 'Hello world!'

process split_letters {
    output:
    path 'chunk_*'

    script:
    """
    printf '${params.str}' | split -b 6 - chunk_
    """
}

process convert_to_upper {
    input:
    path x

    output:
    stdout

    script:
    """
    cat $x | tr '[a-z]' '[A-Z]'
    """
}

workflow {
    split_letters | flatten | convert_to_upper | view { v -> v.trim() }
}
