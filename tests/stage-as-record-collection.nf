#!/usr/bin/env nextflow

nextflow.enable.types = true

record Vcf {
    id: String
    vcf: Path
    tbi: Path?
}

/*
 * files nested in a collection of records should be staged
 * with an explicit stage directive, skipping null values of
 * nullable fields
 */
process combine_vcfs {
    debug true

    input:
    record(
        id: String,
        vcfs: List<Vcf>
    )

    stage:
    stageAs vcfs*.vcf, '*'
    stageAs vcfs*.tbi.findAll { tbi -> tbi != null }, '*'

    script:
    """
    for f in ${vcfs*.vcf.join(' ')} ${vcfs*.tbi.findAll { tbi -> tbi != null }.join(' ')}; do
        [[ -e \$(basename \$f) ]] || { echo "file not staged: \$f"; exit 1; }
    done
    echo ${id}: \$(ls *.vcf *.tbi | sort | xargs)
    """
}

workflow {
    a = file("${workDir}/a.vcf")
    a.text = 'A'
    b = file("${workDir}/b.vcf")
    b.text = 'B'
    b_tbi = file("${workDir}/b.vcf.tbi")
    b_tbi.text = 'B_TBI'

    vcfs = [
        record(id: 'a', vcf: a, tbi: null),
        record(id: 'b', vcf: b, tbi: b_tbi)
    ]
    combine_vcfs( channel.of(record(id: 'all', vcfs: vcfs)) )
}
