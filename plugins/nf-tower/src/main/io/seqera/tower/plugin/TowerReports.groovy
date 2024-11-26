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

import java.nio.charset.Charset
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.PathMatcher
import java.nio.file.Paths
import java.nio.file.StandardCopyOption
import java.time.Duration
import java.util.concurrent.atomic.AtomicInteger

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import groovy.yaml.YamlRuntimeException
import groovy.yaml.YamlSlurper
import groovyx.gpars.agent.Agent
import nextflow.Session
import nextflow.file.FileHelper
import nextflow.trace.GraphObserver
import nextflow.trace.ReportObserver
import nextflow.trace.TimelineObserver
import nextflow.trace.TraceFileObserver
/**
 * If reports are defined at `nf-<workflow_id>-tower.yml`, collects all published files
 * that are reports and writes `nf-<workflow_id>-reports.tsv` file with all the paths.
 *
 * @author Jordi Deu-Pons
 */
@Slf4j
@CompileStatic
class TowerReports {

    private Session session
    private Path launchReportsPath
    private Path workReportsPath
    private PrintWriter reportsFile
    private boolean processReports
    private Agent<PrintWriter> writer
    private List<Map.Entry<String, Map<String, String>>> reportsEntries
    private List<PathMatcher> matchers
    private YamlSlurper yamlSlurper
    private Timer timer
    private AtomicInteger totalReports

    TowerReports(Session session) {
        this.session = session
        this.yamlSlurper = new YamlSlurper()
        this.totalReports = new AtomicInteger(0)
    }

    /**
     * On flow create if there is a tower config yaml for current workflow
     * start writing a reports file at the background.
     *
     * @param workflowId Tower workflow ID
     */
    void flowCreate(String workflowId) {
        Path launchDir = getLaunchDir()
        reportsEntries = parseReportEntries(launchDir, workflowId)
        processReports = reportsEntries.size() > 0
        if (processReports) {
            final fileName = System.getenv().getOrDefault("TOWER_REPORTS_FILE", "nf-${workflowId}-reports.tsv".toString()) as String
            this.launchReportsPath = launchDir.resolve(fileName)
            this.workReportsPath = session?.workDir?.resolve(fileName)
            this.reportsFile = new PrintWriter(Files.newBufferedWriter(launchReportsPath, Charset.defaultCharset()), true)
            this.writer = new Agent<PrintWriter>(reportsFile)

            // send header
            this.writer.send { PrintWriter it -> it.println("key\tpath\tsize\tdisplay\tmime_type") }

            // Schedule a reports copy if launchDir and workDir are different
            if (this.workReportsPath && this.launchReportsPath != this.workReportsPath) {
                final lastTotalReports = new AtomicInteger(0)
                final task = {
                    // Copy the file only if there are new reports
                    if (lastTotalReports.get() < this.totalReports.get()) {
                        try {
                            final total = this.totalReports.get()
                            log.trace("Reports file sync to workdir with ${total} reports")
                            FileHelper.copyPath(launchReportsPath, workReportsPath, StandardCopyOption.REPLACE_EXISTING)
                            lastTotalReports.set(total)
                        } catch (IOException e) {
                            log.error("Error copying reports file ${launchReportsPath.toUriString()} to the workdir ${workReportsPath.toUriString()} -- ${e.message}")
                        }
                    }
                }

                // Copy maximum 1 time per minute
                final oneMinute = Duration.ofMinutes(1).toMillis()
                this.timer = new Timer()
                this.timer.schedule(task, oneMinute, oneMinute)
            }
        }
    }

    /**
     * Retrieve current launch dir
     *
     * @return Launch directory path
     */
    protected Path getLaunchDir() {
        return Paths.get('.').toRealPath()
    }

    /**
     * On flow complete stop writing the reports file at background.
     */
    void flowComplete() {
        if (processReports) {
            if (timer) {
                timer.cancel()
            }
            writer.await()
            // close and upload it
            reportsFile.flush()
            reportsFile.close()
            saveReportsFileUpload()
        }
    }

    protected void saveReportsFileUpload() {
        try {
            FileHelper.copyPath(launchReportsPath, workReportsPath, StandardCopyOption.REPLACE_EXISTING)
            log.debug "Saved reports file ${workReportsPath.toUriString()}"
        }
        catch (Exception e) {
            log.error("Error copying reports file ${launchReportsPath.toUriString()} to the workdir ${workReportsPath.toUriString()} -- ${e.message}")
        }
    }

    /**
     * Load all report entries from tower config yaml file.
     *
     * @param launchDir Nextflow launch directory
     * @param workflowId Tower workflow ID
     */
    protected List<Map.Entry<String, Map<String, String>>> parseReportEntries(Path launchDir, String workflowId) {
        Path towerConfigPath = launchDir.resolve("nf-${workflowId}-tower.yml")

        // Check if Tower config file is define at assets
        if (!Files.exists(towerConfigPath)) {
            towerConfigPath = this.session?.baseDir?.resolve("tower.yml")
        }

        // Load reports definitions if available
        List<Map.Entry<String, Map<String, String>>> reportsEntries = []
        if (towerConfigPath && Files.exists(towerConfigPath)) {
            try {
                final towerConfig = this.yamlSlurper.parse(towerConfigPath)
                if (towerConfig instanceof Map && towerConfig.containsKey("reports")) {
                    Map<String, Map<String, String>> reports = (Map<String, Map<String, String>>) towerConfig.get("reports")
                    for (final e : reports) {
                        reportsEntries.add(e)
                    }
                }
            } catch (YamlRuntimeException e) {
                final msg = e?.cause?.message ?: e.message
                throw new IllegalArgumentException("Invalid tower.yml format -- ${msg}")
            }
        }

        return reportsEntries
    }

    /**
     * On file publish check if the path matches a report pattern a write it to
     * the reports file.
     *
     * @param destination Path of the published file at destination filesystem.
     */
    boolean filePublish(Path destination) {
        if (processReports && destination) {

            if (!matchers) {
                // Initialize report matchers on first event to use the
                // path matcher of the destination filesystem
                matchers = reportsEntries.collect {FileHelper.getPathMatcherFor(convertToGlobPattern(it.key), destination.fileSystem) as PathMatcher}
            }

            for (int p=0; p < matchers.size(); p++) {
                if (matchers.get(p).matches(destination)) {
                    final reportEntry = this.reportsEntries.get(p)
                    writer.send((PrintWriter it) -> writeRecord(it, reportEntry, destination))
                    return true
                }
            }
        }
        return false
    }

    private writeRecord(PrintWriter it, Map.Entry<String,Map<String,String>> reportEntry, Path destination) {
        try {
            final target = destination.toUriString()
            final numRep = totalReports.incrementAndGet()
            log.trace("Adding report [${numRep}] ${reportEntry.key} -- ${target}")
            // Report properties
            final display = reportEntry.value.get("display", "")
            final mimeType = reportEntry.value.get("mimeType", "")
            it.println("${reportEntry.key}\t${target}\t${destination.size()}\t${display}\t${mimeType}");
            it.flush()
        }
        catch (Throwable e) {
            log.error ("Unexpected error writing report entry '${destination.toUriString()}'", e)
        }
    }

    protected static String convertToGlobPattern(String reportKey) {
        final prefix = reportKey.startsWith("**/") ? "" : "**/"
        return "glob:${prefix}${reportKey}"
    }

    /**
     * Publish any built-in reports that are enabled to the reports file.
     */
    void publishRuntimeReports() {
        final config = session.config
        final files = []

        if( config.navigate('report.enabled') )
            files << config.navigate('report.file', ReportObserver.DEF_FILE_NAME)

        if( config.navigate('timeline.enabled') )
            files << config.navigate('timeline.file', TimelineObserver.DEF_FILE_NAME)

        if( config.navigate('trace.enabled') )
            files << config.navigate('trace.file', TraceFileObserver.DEF_FILE_NAME)

        if( config.navigate('dag.enabled') )
            files << config.navigate('dag.file', GraphObserver.DEF_FILE_NAME)

        for( def file : files )
            filePublish( (file as Path).complete() )
    }
}
