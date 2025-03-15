process python_version {
    array 3
    input:
    val x
    path input_file

    output:
    stdout

    script:
    """
    echo $x
    cat $input_file
    """
}

workflow {
    python_version(Channel.of(1,2,3), params.array_file_path) | view
}

