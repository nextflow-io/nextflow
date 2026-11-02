/*
 * End-to-end test of Nextflow running natively on Windows, with every
 * task executed in a Linux container. See validation/windows/test.ps1
 */
params.input = "$projectDir/data/greeting.txt"
params.notes = "$projectDir/data/my data/notes.txt"
params.image = 'public.ecr.aws/docker/library/ubuntu:24.04'

// input file and script from the project `bin` directory
process UPPER {
    container params.image
    publishDir 'results', mode: 'copy'

    input:
    path x

    output:
    path 'upper.txt'

    script:
    """
    shout.sh $x > upper.txt
    """
}

// input file in a path with blanks, scratch directory, env output
// and publishing with the default (symlink) mode
process COUNT {
    container params.image
    scratch true
    publishDir 'results'

    input:
    path x
    path notes

    output:
    path 'count.txt'
    env 'LINES'

    script:
    """
    cat $x $notes | wc -l > count.txt
    LINES=\$(wc -l < $x)
    """
}

// report where the task actually runs
process ENV_INFO {
    container params.image

    output:
    stdout

    script:
    """
    echo "kernel: \$(uname -s); workdir: \$PWD"
    """
}

// failed task retried
process RETRY {
    container params.image
    errorStrategy 'retry'
    maxRetries 1

    input:
    val n

    output:
    stdout

    script:
    """
    if [ $task.attempt -lt 2 ]; then echo 'failing on purpose' >&2; exit 3; fi
    echo "attempt $task.attempt ok"
    """
}

// task container stopped when the time limit is exceeded
process SLOW {
    container params.image
    time '10s'
    errorStrategy 'ignore'

    script:
    """
    sleep 300
    """
}

workflow {
    upper = UPPER(file(params.input))
    COUNT(upper, file(params.notes))
    COUNT.out[1].view { n -> "lines: ${n.trim()}" }
    ENV_INFO().view { s -> s.trim() }
    RETRY(1).view { s -> s.trim() }
    SLOW()
}
