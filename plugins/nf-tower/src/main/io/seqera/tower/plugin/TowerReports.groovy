package io.seqera.tower.plugin

import groovy.yaml.YamlSlurper
import groovyx.gpars.agent.Agent
import nextflow.file.FileHelper

import java.nio.charset.Charset
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.PathMatcher
import java.nio.file.Paths
import java.util.stream.Collectors


/**
 * If reports are defined at `nf-<workflow_id>-tower.yml`, collects all published files
 * that are reports and writes `nf-<workflow_id>-reports.tsv` file with all the paths.
 *
 * @author Jordi Deu-Pons
 */
class TowerReports {

    private Path reportsPath
    private PrintWriter reportsFile
    private Agent<PrintWriter> reportsWriter
    private List<String> reportPatterns
    private List<PathMatcher> reportMatchers
    private YamlSlurper yamlSlurper

    TowerReports() {
        this.yamlSlurper = new YamlSlurper()
    }

    /**
     * On flow create if there is a tower config yaml for current workflow
     * start writing a reports file at the background.
     *
     * @param workflowId Tower workflow ID
     */
    void flowCreate(String workflowId) {
        Path launchDir = Paths.get('.').toRealPath()
        loadReportPatterns(launchDir, workflowId)
        if (reportPatterns) {
            this.reportsPath = launchDir.resolve(System.getenv().getOrDefault("TOWER_REPORTS_FILE", "nf-${workflowId}-reports.tsv"))
            this.reportsFile = new PrintWriter(Files.newBufferedWriter(reportsPath, Charset.defaultCharset()))
            this.reportsWriter = new Agent<PrintWriter>(reportsFile)

            // send header
            this.reportsWriter.send { PrintWriter it -> it.println("key\tpath")}
        }
    }

    /**
     * On flow complete stop writing the reports file at background.
     */
    void flowComplete() {
        if (reportsWriter) {
            reportsWriter.await()
            reportsFile.flush()
            reportsFile.close()
        }
    }

    /**
     * Load all report glob patterns from tower config yaml file.
     *
     * @param launchDir Nextflow launch directory
     * @param workflowId Tower workflow ID
     */
    protected void loadReportPatterns(Path launchDir, String workflowId) {
        Path towerConfigPath = launchDir.resolve("nf-${workflowId}-tower.yml")
        if (Files.exists(towerConfigPath)) {
            final towerConfig = yamlSlurper.parse(towerConfigPath)
            this.reportPatterns = new ArrayList<>()
            if (towerConfig instanceof Map && towerConfig.containsKey("reports")) {
                final reports = towerConfig.get("reports") as Map
                reports.keySet().forEach(k -> {
                    if (k) {
                        reportPatterns.add(k.toString())
                    }
                })
            }
        }
    }

    /**
     * On file publish check if the path matches a report pattern a write it to
     * the reports file.
     *
     * @param destination Path of the published file at destination filesystem.
     */
    void filePublish(Path destination) {
        if (reportsWriter && destination) {

            if (!reportMatchers) {
                // Initialize report matchers on first event to use the
                // path matcher of the destination filesystem
                reportMatchers = reportPatterns.stream()
                        .map(p -> FileHelper.getPathMatcherFor("glob:**/${p}", destination.fileSystem))
                        .collect(Collectors.toList())
            }

            for (int p=0; p < reportMatchers.size(); p++) {
                if (reportMatchers.get(p).matches(destination)) {
                    final dst = destination.toAbsolutePath().toString()
                    final pattern = reportPatterns.get(p)
                    reportsWriter.send { PrintWriter it -> it.println("${pattern}\t${dst}") }
                    break
                }
            }
        }
    }
}
