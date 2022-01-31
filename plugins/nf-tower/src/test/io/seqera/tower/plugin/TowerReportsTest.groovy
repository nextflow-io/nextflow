package io.seqera.tower.plugin

import nextflow.Session
import spock.lang.Specification

/**
 *
 * @author Jordi Deu-Pons <jordi@seqera.io>
 */
class TowerReportsTest extends Specification {

    def 'should convert to glob pattern'() {

        expect:
        TowerReports.convertToGlobPattern(KEY) == EXPECTED

        where:
        KEY               | EXPECTED
        "multiqc.html"    | "glob:**/multiqc.html"
        "**/multiqc.html" | "glob:**/multiqc.html"
        "*/multiqc.html"  | "glob:**/*/multiqc.html"
        "reports/*.html"  | "glob:**/reports/*.html"

    }

    def 'should load reports from tower yml file'() {

        given:
        def launchDir = File.createTempDir()
        def workflowId = "1khBUM1SUioskd"
        def config = new File(launchDir, "nf-${workflowId}-tower.yml")
        config.text = """
        reports:
          multiqc_report.html:
            display: "MultiQC HTML report"
          deseq2.plots.pdf:
            display: "All samples STAR Salmon DESeq2 QC PDF plots"
          salmon.merged.gene_counts.tsv:
            display: "All samples STAR Salmon merged gene raw counts"
          "*.merged.gene_counts.tsv":
            display: "All samples STAR Salmon merged gene raw counts"
        """.stripIndent()

        when:
        def reports = new TowerReports(Mock(Session))
        def entries = reports.parseReportEntries(launchDir.toPath(), workflowId)

        then:
        entries == [
                "multiqc_report.html"          : ["display": "MultiQC HTML report"],
                "deseq2.plots.pdf"             : ["display": "All samples STAR Salmon DESeq2 QC PDF plots"],
                "salmon.merged.gene_counts.tsv": ["display": "All samples STAR Salmon merged gene raw counts"],
                "*.merged.gene_counts.tsv"     : ["display": "All samples STAR Salmon merged gene raw counts"]
        ].collect()


    }


}
