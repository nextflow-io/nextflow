/*
 * Copyright 2013-2024, Seqera Labs
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

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

    def 'should generate reports file'() {

        given: 'a launch directory with a tower yaml file'
        def launchDir = File.createTempDir()
        def workflowId = "1khBUM1SUioskd"
        def config = new File(launchDir, "nf-${workflowId}-tower.yml")
        config.text = """
        reports:
          multiqc_report.html:
            display: "MultiQC HTML report"
          deseq2.plots.pdf:
            display: "All samples STAR Salmon DESeq2 QC PDF plots"
            mimeType: "application/pdf"
          "*.merged.gene_counts.tsv":
            display: "All samples STAR Salmon merged gene raw counts"
        """.stripIndent()

        and: 'a tower reports instance'
        def session = new Session()
        TowerReports reports = Spy(TowerReports, constructorArgs: [session])
        reports.launchDir >> launchDir.toPath()

        and: 'some reports'
        def repo1 = new File(launchDir, "multiqc_report.html")
        repo1.text = "html"
        def repo2 = new File(launchDir, "deseq2.plots.pdf")
        repo2.text = "pdf"
        def repo3 = new File(launchDir, "salmon.merged.gene_counts.tsv")
        repo3.text = "tsv"

        when: 'a workflow runs'
        reports.flowCreate(workflowId)
        reports.filePublish(repo1.toPath())
        reports.filePublish(repo2.toPath())
        reports.filePublish(repo3.toPath())
        reports.flowComplete()

        then:
        def result = new File(launchDir, "nf-${workflowId}-reports.tsv")
        result.text == "key\tpath\tsize\tdisplay\tmime_type\n" +
                "multiqc_report.html\t${repo1.toPath().toUriString()}\t4\tMultiQC HTML report\t\n" +
                "deseq2.plots.pdf\t${repo2.toPath().toUriString()}\t3\tAll samples STAR Salmon DESeq2 QC PDF plots\tapplication/pdf\n" +
                "*.merged.gene_counts.tsv\t${repo3.toPath().toUriString()}\t3\tAll samples STAR Salmon merged gene raw counts\t\n"

    }


}
