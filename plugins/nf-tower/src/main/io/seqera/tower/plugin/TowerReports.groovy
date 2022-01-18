package io.seqera.tower.plugin

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import groovy.yaml.YamlSlurper
import groovyx.gpars.agent.Agent
import nextflow.Global
import nextflow.file.FileHelper

import java.nio.charset.Charset
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.PathMatcher
import java.nio.file.Paths
import java.time.Duration
import java.util.concurrent.atomic.AtomicInteger
import java.util.stream.Collectors

/**
 * If reports are defined at `nf-<workflow_id>-tower.yml`, collects all published files
 * that are reports and writes `nf-<workflow_id>-reports.tsv` file with all the paths.
 *
 * @author Jordi Deu-Pons
 */
@Slf4j
@CompileStatic
class TowerReports {

    private Path launchReportsPath
    private Path workReportsPath
    private PrintWriter reportsFile
    private boolean processReports
    private Agent<PrintWriter> writer
    private List<String> patterns
    private List<PathMatcher> matchers
    private YamlSlurper yamlSlurper
    private Timer timer
    private AtomicInteger totalReports

    TowerReports() {
        this.yamlSlurper = new YamlSlurper()
        this.timer = new Timer()
        this.totalReports = new AtomicInteger(0)
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
        if (processReports) {
            final fileName = System.getenv().getOrDefault("TOWER_REPORTS_FILE", "nf-${workflowId}-reports.tsv")
            this.launchReportsPath = launchDir.resolve(fileName)
            this.workReportsPath = Global?.session?.workDir?.resolve(fileName)
            this.reportsFile = new PrintWriter(Files.newBufferedWriter(launchReportsPath, Charset.defaultCharset()), true)
            this.writer = new Agent<PrintWriter>(reportsFile)

            // send header
            this.writer.send { PrintWriter it -> it.println("key\tpath\tsize")}

            // Schedule a reports copy if launchDir and workDir are different
            if (this.workReportsPath && this.launchReportsPath != this.workReportsPath) {
                final lastTotalReports = new AtomicInteger(0)
                final task = {
                    // Copy the file only if there are new reports
                    if (lastTotalReports.get() < this.totalReports.get()) {
                        try {
                            final total = this.totalReports.get()
                            FileHelper.copyPath(launchReportsPath, workReportsPath)
                            lastTotalReports.set(total)
                        } catch (IOException e) {
                            log.error("Copying reports file {} to the workdir.", this.launchReportsPath)
                        }
                    }
                }

                // Copy maximum 1 time per minute
                final oneMinute = Duration.ofMinutes(1).toMillis()
                this.timer.schedule(task, oneMinute, oneMinute)
            }
        }
    }

    /**
     * On flow complete stop writing the reports file at background.
     */
    void flowComplete() {
        if (processReports) {
            timer.cancel()
            writer.await()
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
        processReports = false
        Path towerConfigPath = launchDir.resolve("nf-${workflowId}-tower.yml")
        if (Files.exists(towerConfigPath)) {
            final towerConfig = yamlSlurper.parse(towerConfigPath)
            this.patterns = new ArrayList<>()
            if (towerConfig instanceof Map && towerConfig.containsKey("reports")) {
                final reports = towerConfig.get("reports") as Map
                reports.keySet().forEach(k -> {
                    if (k) {
                        patterns.add(k.toString())
                    }
                })
            }
            processReports = this.patterns.size() > 0
        }
    }

    /**
     * On file publish check if the path matches a report pattern a write it to
     * the reports file.
     *
     * @param destination Path of the published file at destination filesystem.
     */
    void filePublish(Path destination) {
        if (processReports && destination) {

            if (!matchers) {
                // Initialize report matchers on first event to use the
                // path matcher of the destination filesystem
                matchers = patterns.stream()
                        .map(p -> FileHelper.getPathMatcherFor("glob:**/${p}", destination.fileSystem))
                        .collect(Collectors.toList())
            }

            for (int p=0; p < matchers.size(); p++) {
                if (matchers.get(p).matches(destination)) {
                    final dst = destination.toUriString()
                    final pattern = patterns.get(p)
                    writer.send { PrintWriter it -> it.println("${pattern}\t${dst}\t${destination.size()}") }
                    final numRep = totalReports.incrementAndGet()
                    log.debug("Adding report [{}] {} / {}", numRep, pattern, dst)
                    break
                }
            }
        }
    }
}
