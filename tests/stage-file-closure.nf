
nextflow.preview.types = true

process make_files {
    input:
    n_groups: Integer
    group_size: Integer

    output:
    files('*/*.txt')

    script:
    """
    for i in `seq 1 ${n_groups}`; do
        mkdir group\${i}
        for j in `seq 1 ${group_size}`; do
            touch group\${i}/sample\${j}.txt
        done
    done
    """
}

process ls {
    input:
    slice: Set<Path>

    stage:
    stageAs(slice) { file -> "${file.parent.name}/${file.name}" }

    output:
    stdout()

    script:
    """
    ls -1 */*.txt | sort
    """
}

workflow {
    ls(make_files(3, 3)).view()
}
