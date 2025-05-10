params.str = 'Hello world!'

process splitLetters {
    output:
    path 'chunk_*'

    script:
    """
    printf '${params.str}' | split -b 6 - chunk_
    """
}

process convertToUpper {
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
    splitLetters | flatten | convertToUpper | view { v -> v.trim() }
}
