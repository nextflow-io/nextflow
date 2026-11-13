nextflow.enable.types = true

process MULTIQC {
    input:
    counts: Bag<Path>

    output:
    file('multiqc_report.html')

    script:
    """
    echo "summarized ${counts.size()} samples" > multiqc_report.html
    """
}
