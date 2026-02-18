nextflow.preview.recursion = true

params.input = "recurse-workflow.in"

workflow {
    clock
        .recurse(file(params.input))
        .until { file -> file.size() > 64 }
        .view { file -> file.text }
}

workflow clock {
    take:
    logfile

    emit:
    tock(tick(logfile))
}

process tick {
    input:
    path 'input.txt'

    output:
    path 'result.txt'

    script:
    """
    cat input.txt > result.txt
    echo "Task ${task.index} : tick" >> result.txt
    """
}

process tock {
    input:
    path 'input.txt'

    output:
    path 'result.txt'

    script:
    """
    cat input.txt > result.txt
    echo "Task ${task.index} : tock" >> result.txt
    """
}
