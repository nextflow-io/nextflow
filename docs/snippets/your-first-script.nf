nextflow.enable.types = true

// Pipeline parameters
params {
    str: String = "Hello world!"
}

// Entry workflow
workflow {
    main:
    ch_str = channel.of(params.str)                         // Create a channel of the pipeline input
    ch_chunks = split(ch_str).flatMap { chunks -> chunks }  // Split string into chunks and emit each chunk separately
    ch_upper = convert_to_upper(ch_chunks)                  // Convert lowercase letters to uppercase letters

    publish:
    lower = ch_chunks
    upper = ch_upper
}

// Pipeline outputs
output {
    lower {
        path 'lower'
    }
    upper {
        path 'upper'
    }
}

// split process
process split {
    input:
    x: String

    output:
    files('chunk_*')

    script:
    """
    printf '${x}' | split -b 6 - chunk_
    """
}

// convert_to_upper process
process convert_to_upper {
    tag "$y"

    input:
    y: Path

    output:
    file("upper_${y}")

    script:
    """
    cat $y | tr '[a-z]' '[A-Z]' > upper_${y}
    """
}
