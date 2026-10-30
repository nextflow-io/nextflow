nextflow.enable.types = true

process MULTIQC {
    input:
    logs: Bag<Path>
    args: String?

    output:
    file('multiqc_report.html')

    script:
    """
    echo "multiqc ${args ?: ''} ${logs.toSorted().join(' ')}" > multiqc_report.html
    """
}
