#!/usr/bin/env nextflow

nextflow.enable.types = true

record Vcf {
    id: String
    vcf: Path
    tbi: Path?
}

/*
 * files nested in a collection of records should be staged
 * implicitly, skipping null values of nullable fields
 */
process combine_vcfs {
    debug true

    input:
    record(
        id: String,
        vcfs: List<Vcf>
    )

    script:
    """
    for f in ${vcfs*.vcf.join(' ')} ${vcfs*.tbi.findAll { tbi -> tbi != null }.join(' ')}; do
        [[ "\$f" == "\$(basename \$f)" ]] || { echo "not a staged path: \$f"; exit 1; }
        [[ -e \$f ]] || { echo "file not staged: \$f"; exit 1; }
    done
    echo ${id}: \$(ls *.vcf *.tbi | sort | xargs)
    """
}

/*
 * a nullable record should allow null values
 */
process index_vcf_opt {
    debug true

    input:
    vcf: Vcf?

    script:
    """
    echo vcf: ${vcf == null ? 'none' : vcf.vcf}
    """
}

/*
 * a nullable record collection should allow null values
 */
process combine_vcfs_opt {
    debug true

    input:
    vcfs: List<Vcf>?

    script:
    """
    echo vcfs: ${vcfs == null ? 'none' : vcfs*.vcf.join(' ')}
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
    combine_vcfs_opt( channel.of(null) )
    index_vcf_opt( channel.of(null) )
}
