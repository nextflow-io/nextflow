nextflow.preview.recursion = true

params.start = 10

workflow {
    count_down
        .recurse(params.start)
        .until { v -> v == 0 }
        .view { v -> "${v}..." }
}

process count_down {
    input:
    val v

    output:
    val v

    exec:
    sleep(1000)
    v = v - 1
}
