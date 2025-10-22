/*
 * Copyright 2013-2025, Seqera Labs
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
 */

package io.seqera.tower.plugin.launch

import groovy.json.JsonBuilder
import groovy.json.JsonSlurper
import groovy.transform.Canonical
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import io.seqera.tower.plugin.BaseCommandImpl
import io.seqera.tower.plugin.TowerClient
import nextflow.BuildInfo
import nextflow.cli.CmdLaunch
import nextflow.util.ColorUtil
import nextflow.exception.AbortOperationException
import nextflow.file.FileHelper
import nextflow.scm.AssetManager
import nextflow.util.SpinnerUtil
import org.yaml.snakeyaml.Yaml

import java.net.http.HttpRequest
import java.net.http.HttpResponse
import java.nio.file.NoSuchFileException
import java.nio.file.Path
import java.util.concurrent.atomic.AtomicBoolean
import java.util.regex.Pattern

/**
 * CLI sub-command LAUNCH -- Launch a workflow in Seqera Platform
 *
 * @author Phil Ewels <phil.ewels@seqera.io>
 */
@Slf4j
@CompileStatic
class LaunchCommandImpl extends BaseCommandImpl implements CmdLaunch.LaunchCommand {

    // ===== Constants =====

    static final public List<String> VALID_PARAMS_FILE = ['json', 'yml', 'yaml']
    static final private Pattern DOT_ESCAPED = ~/\\\./
    static final private Pattern DOT_NOT_ESCAPED = ~/(?<!\\)\./
    static final int LOG_POLL_INTERVAL_MS = 2_000
    static final int LOG_GRACE_PERIOD_MS = 5_000
    static final int LOG_SLEEP_INTERVAL_MS = 100

    // ===== Data Classes =====

    /**
     * Holds all context needed for launching a workflow
     */
    @Canonical
    static class LaunchContext {
        String accessToken
        String apiEndpoint
        String userName
        Long workspaceId
        String orgName
        String workspaceName
        String computeEnvId
        String computeEnvName
        String workDir
    }

    /**
     * Holds the result of a workflow launch
     */
    @Canonical
    static class WorkflowLaunchResult {
        String workflowId
        String runName
        String commitId
        String revision
        String repository
        String trackingUrl
    }

    // ===== Main Entry Point =====

    @Override
    void launch(CmdLaunch.LaunchOptions options) {
        printBanner(options)

        // Validate and resolve pipeline
        final resolvedPipelineUrl = validateAndResolvePipeline(options.pipeline)

        // Initialize launch context with auth and workspace info
        final context = initializeLaunchContext(options)

        // Build and submit the launch request
        final result = submitWorkflowLaunch(options, context, resolvedPipelineUrl)

        // Display launch information
        printLaunchInfo(result.repository, result.runName, result.commitId, result.revision,
                        context.workDir, context.computeEnvName, context.userName,
                        context.orgName, context.workspaceName, options)
        printSuccessMessage(result.workflowId, result.trackingUrl, options)

        // Poll for workflow logs
        if (result.workflowId) {
            pollWorkflowLogs(result.workflowId, context.workspaceId, result.trackingUrl,
                           context.accessToken, context.apiEndpoint, options)
        }
    }

    // ===== Launch Phases =====

    /**
     * Validate pipeline path and resolve to full repository URL
     */
    private String validateAndResolvePipeline(String pipeline) {
        log.debug "Pipeline repository: ${pipeline}"

        if (isLocalPath(pipeline)) {
            log.debug "Rejecting local file path: ${pipeline}"
            throw new AbortOperationException("Local file paths are not supported. Please provide a remote repository URL.")
        }

        final resolvedUrl = resolvePipelineUrl(pipeline)
        log.debug "Resolved pipeline URL: ${resolvedUrl}"
        return resolvedUrl
    }

    /**
     * Initialize launch context by loading config and resolving workspace/compute environment
     */
    private LaunchContext initializeLaunchContext(CmdLaunch.LaunchOptions options) {
        log.debug "Initializing launch context"

        // Load configuration
        final config = readConfig()
        final apiEndpoint = (config['tower.endpoint'] ?: TowerClient.DEF_ENDPOINT_URL) as String
        final accessToken = config['tower.accessToken'] as String

        if (!accessToken) {
            throw new AbortOperationException("No authentication found. Please run 'nextflow auth login' first.")
        }

        // Resolve workspace
        final workspaceId = resolveWorkspaceId(config, options.workspace, accessToken, apiEndpoint)
        final userName = getUserInfo(accessToken, apiEndpoint).name as String

        String orgName = null
        String workspaceName = null
        if (workspaceId) {
            final wsDetails = getWorkspaceDetails(accessToken, apiEndpoint, workspaceId.toString())
            orgName = wsDetails?.orgName as String
            workspaceName = wsDetails?.workspaceName as String
            log.debug "Using workspace '${workspaceName}' (ID: ${workspaceId})"
        } else {
            log.debug "Using personal workspace for user: ${userName}"
        }

        // Resolve compute environment
        final computeEnvInfo = resolveComputeEnvironment(config, options.computeEnv, workspaceId, accessToken, apiEndpoint)
        final workDir = resolveWorkDirectory(options.workDir, computeEnvInfo)

        return new LaunchContext(
            accessToken: accessToken,
            apiEndpoint: apiEndpoint,
            userName: userName,
            workspaceId: workspaceId,
            orgName: orgName,
            workspaceName: workspaceName,
            computeEnvId: computeEnvInfo.id as String,
            computeEnvName: computeEnvInfo.name as String,
            workDir: workDir
        )
    }

    /**
     * Build launch request, submit to API, and return result
     */
    private WorkflowLaunchResult submitWorkflowLaunch(CmdLaunch.LaunchOptions options, LaunchContext context, String pipelineUrl) {
        log.debug "Submitting workflow launch"

        // Build request payload
        final paramsText = buildParamsText(options.params, options.paramsFile)
        final configText = buildConfigText(options.configFiles)
        final launchRequest = buildLaunchRequestPayload(options, context, pipelineUrl, paramsText, configText)

        // Submit to API
        final queryParams = context.workspaceId ? [workspaceId: context.workspaceId.toString()] : [:]
        final response = apiPost('/workflow/launch', launchRequest, queryParams, context.accessToken, context.apiEndpoint)

        // Fetch workflow details for accurate launch info
        final workflowDetails = fetchWorkflowDetails(response.workflowId as String, context.workspaceId,
                                                      context.accessToken, context.apiEndpoint)

        // Extract and return launch result
        return extractLaunchResult(response, workflowDetails, options, pipelineUrl, context)
    }

    /**
     * Build the launch request payload
     */
    private Map buildLaunchRequestPayload(CmdLaunch.LaunchOptions options, LaunchContext context,
                                          String pipelineUrl, String paramsText, String configText) {
        def launch = [:]
        launch.computeEnvId = context.computeEnvId
        launch.workDir = context.workDir
        launch.pipeline = pipelineUrl
        launch.resume = options.resume != null
        launch.pullLatest = options.latest
        launch.stubRun = options.stubRun

        if (options.runName) launch.runName = options.runName
        if (options.revision) launch.revision = options.revision
        if (options.profile) launch.configProfiles = options.profile
        if (configText) launch.configText = configText
        if (paramsText) launch.paramsText = paramsText
        if (options.mainScript) launch.mainScript = options.mainScript
        if (options.entryName) launch.entryName = options.entryName

        log.debug "Built launch request with ${launch.size()} parameters"
        return [launch: launch]
    }

    /**
     * Fetch workflow details from API if workflow ID is available
     */
    private Map fetchWorkflowDetails(String workflowId, Long workspaceId, String accessToken, String apiEndpoint) {
        if (!workflowId) return null

        log.debug "Fetching workflow details for ID: ${workflowId}"
        final queryParams = workspaceId ? [workspaceId: workspaceId.toString()] : [:]
        return apiGet("/workflow/${workflowId}", queryParams, accessToken, apiEndpoint)
    }

    /**
     * Extract launch result from API response and workflow details
     */
    private WorkflowLaunchResult extractLaunchResult(Map response, Map workflowDetails,
                                                     CmdLaunch.LaunchOptions options, String pipelineUrl, LaunchContext context) {
        def runName = 'unknown'
        def commitId = 'unknown'
        def revision = options.revision

        if (workflowDetails?.workflow) {
            final workflow = workflowDetails.workflow as Map
            runName = workflow.runName as String ?: options.runName ?: 'unknown'
            commitId = workflow.commitId as String ?: 'unknown'
            revision = workflow.revision as String ?: options.revision
        } else {
            runName = options.runName ?: 'unknown'
        }

        // Build tracking URL
        final webUrl = getWebUrlFromApiEndpoint(context.apiEndpoint)
        String trackingUrl = null
        if (response.workflowId) {
            if (context.orgName && context.workspaceName) {
                trackingUrl = "${webUrl}/orgs/${context.orgName}/workspaces/${context.workspaceName}/watch/${response.workflowId}/"
            } else if (context.userName) {
                trackingUrl = "${webUrl}/user/${context.userName}/watch/${response.workflowId}/"
            }
        }

        return new WorkflowLaunchResult(
            workflowId: response.workflowId as String,
            runName: runName,
            commitId: commitId,
            revision: revision,
            repository: pipelineUrl,
            trackingUrl: trackingUrl
        )
    }

    // ===== Configuration & Resolution Methods =====

    /**
     * Resolve compute environment by flag name config computeEnvId or get primary
     */
    protected Map resolveComputeEnvironment(Map config, String computeEnvName, Long workspaceId, String accessToken, String apiEndpoint) {
        Map computeEnvInfo = null
        if (!computeEnvName && config?.get('tower.computeEnvId')) {
            computeEnvInfo = getComputeEnvironment(accessToken, apiEndpoint, config['tower.computeEnvId'] as String, workspaceId?.toString())
        } else {
            log.debug "Looking up compute environment: ${computeEnvName ?: '(primary)'}"
            computeEnvInfo = findComputeEnv(computeEnvName, workspaceId, accessToken, apiEndpoint)
        }
        if (!computeEnvInfo) {
            if (computeEnvName) {
                throw new AbortOperationException("Compute environment '${computeEnvName}' not found")
            } else {
                throw new AbortOperationException("No primary compute environment found")
            }
        }

        log.debug "Using compute environment: '${computeEnvInfo.name}' (ID: ${computeEnvInfo.id})"
        return computeEnvInfo
    }

    /**
     * Resolve work directory from CLI option or compute environment
     */
    private String resolveWorkDirectory(String cliWorkDir, Map computeEnvInfo) {
        String workDir = cliWorkDir
        if (!workDir && computeEnvInfo.workDir) {
            workDir = computeEnvInfo.workDir as String
        }
        if (!workDir) {
            throw new AbortOperationException("Work directory is required. Please specify -w/--work-dir or ensure your compute environment has a workDir configured.")
        }
        return workDir
    }

    // ===== Display Methods =====

    protected void printBanner(CmdLaunch.LaunchOptions options) {
        if (ColorUtil.isAnsiEnabled()) {
            // Plain header for verbose log
            log.debug "N E X T F L O W  ~  version ${BuildInfo.version}"

            // Fancy coloured header for the console output
            println ""
            // Use exact colour codes matching the Nextflow green brand
            final BACKGROUND = "\033[1m\033[38;5;232m\033[48;5;43m"
            print("$BACKGROUND N E X T F L O W ")
            print("\033[0m")

            // Show Nextflow version and launch context
            print(ColorUtil.colorize("  ~  ", "dim", true))
            println("version " + BuildInfo.version)
            println("Launching workflow in Seqera Platform")
            println ""
        } else {
            // Plain header to the console if ANSI is disabled
            log.info "N E X T F L O W  ~  version ${BuildInfo.version}"
            log.info "Launching workflow in Seqera Platform"
        }
    }

    protected void printLaunchInfo(String repo, String runName, String commitId, String revision, String workDir, String computeEnvName, String userName, String orgName, String workspaceName, CmdLaunch.LaunchOptions options) {
        def showRevision = commitId && commitId != 'unknown'
        def showRevisionBrackets = revision && revision != 'unknown'

        if (ColorUtil.isAnsiEnabled()) {
            def debugMsg = "Launched `${repo}` [${runName}]"
            if (showRevision) {
                debugMsg += " - revision: ${commitId}"
            }
            if (showRevisionBrackets) {
                debugMsg += " [${revision}]"
            }
            log.debug debugMsg

            print("Launched")
            print(ColorUtil.colorize(" `$repo` ", "magenta", true))
            print(ColorUtil.colorize("[", "dim", true))
            print(ColorUtil.colorize(runName, "cyan bold", true))
            print(ColorUtil.colorize("]", "dim", true))
            if (showRevision) {
                print(" - ")
                print(ColorUtil.colorize("revision: $commitId", "cyan", true))
            }
            if (showRevisionBrackets) {
                print(ColorUtil.colorize(" [", "dim", true))
                print(ColorUtil.colorize(revision, "cyan", true))
                print(ColorUtil.colorize("]", "dim", true))
            }
            println ""

            // Print username
            if (userName) {
                println(" 👤 user: ${ColorUtil.colorize(userName, 'cyan', true)}")
            }

            // Print workspace
            if (orgName && workspaceName) {
                println(" 🏢 workspace: ${ColorUtil.colorize(orgName + ' / ' + workspaceName, 'cyan', true)}")
            } else {
                println(" 🏢 workspace: ${ColorUtil.colorize("Personal workspace", 'cyan', true)}")
            }

            // Print work directory
            println(" 📁 workdir: ${ColorUtil.colorize(workDir, 'cyan', true)}")

            // Print compute environment
            println(" ☁️ compute: ${ColorUtil.colorize(computeEnvName, 'cyan', true)}\n")
        } else {
            def plainMsg = "Launched `${repo}` [${runName}]"
            if (showRevision) {
                plainMsg += " - revision: ${commitId}"
            }
            if (showRevisionBrackets) {
                plainMsg += " [${revision}]"
            }
            log.info plainMsg
            if (userName) {
                log.info " user: ${userName}"
            }
            if (orgName && workspaceName) {
                log.info " workspace: ${orgName} / ${workspaceName}"
            } else {
                log.info " workspace: Personal workspace"
            }
            log.info " workdir: ${workDir}"
            log.info " compute: ${computeEnvName}"
        }
    }

    protected void printSuccessMessage(String workflowId, String trackingUrl, CmdLaunch.LaunchOptions options) {
        if (ColorUtil.isAnsiEnabled()) {
            if (trackingUrl) {
                print(ColorUtil.colorize("Workflow launched successfully: ", "green"))
                ColorUtil.printColored(trackingUrl, "cyan")
            } else if (workflowId) {
                ColorUtil.printColored("Workflow launched successfully!", "green")
                print("Workflow ID: ")
                ColorUtil.printColored(workflowId, "cyan")
            } else {
                ColorUtil.printColored("Workflow launched successfully!", "green")
            }
            print(ColorUtil.colorize("💡 To do more with Seqera Platform via the CLI, see ", "dim"))
            ColorUtil.printColored("https://github.com/seqeralabs/tower-cli/", "cyan bold")
            println ""
        } else {
            if (trackingUrl) {
                log.info "Workflow launched successfully: ${trackingUrl}"
            } else if (workflowId) {
                log.info "Workflow launched successfully!"
                log.info "Workflow ID: ${workflowId}"
            } else {
                log.info "Workflow launched successfully!"
            }
            log.info "To do more with Seqera Platform via the CLI, see https://github.com/seqeralabs/tower-cli/"
        }
    }

    // ===== Log Polling Methods =====

    /**
     * Poll workflow logs until the workflow completes
     */
    private void pollWorkflowLogs(String workflowId, Long workspaceId, String trackingUrl,
                                  String accessToken, String apiEndpoint, CmdLaunch.LaunchOptions options) {
        log.debug "Starting log polling for workflow ID: ${workflowId}"

        final queryParams = workspaceId ? [workspaceId: workspaceId.toString()] : [:]
        final shouldExit = new AtomicBoolean(false)
        final spinner = new SpinnerUtil("Checking workflow status...", true) // Use waiting animation initially
        final shutdownHook = createShutdownHook(shouldExit, trackingUrl, spinner)

        Runtime.getRuntime().addShutdownHook(shutdownHook)

        log.info "Workflow submitted, awaiting log output. It is now safe to exit with ctrl+c\n"
        println "\n" + ('─' * 70)

        try {
            runLogPollingLoop(workflowId, queryParams, accessToken, apiEndpoint, shouldExit, spinner)
        } finally {
            removeShutdownHook(shutdownHook)
        }

        log.debug "Log polling completed for workflow ID: ${workflowId}"
    }

    /**
     * Create shutdown hook for graceful exit on Ctrl+C
     */
    private Thread createShutdownHook(AtomicBoolean shouldExit, String trackingUrl, SpinnerUtil spinner) {
        return new Thread({
            shouldExit.set(true)
            log.debug "Shutdown hook triggered - setting exit flag"

            // Stop spinner if running
            if (spinner.isRunning()) {
                Thread.sleep(50) // Give spinner thread a moment to stop
                spinner.stop()
            }

            printLogPollingExitMessage(trackingUrl)
        })
    }

    /**
     * Run the main log polling loop
     */
    private void runLogPollingLoop(String workflowId, Map queryParams, String accessToken,
                                    String apiEndpoint, AtomicBoolean shouldExit, SpinnerUtil spinner) {
        int displayedLogCount = 0
        def firstLogReceived = false
        def workflowFinished = false
        def workflowFinishedTime = 0L
        final finalStatuses = Set.of('SUCCEEDED', 'FAILED', 'CANCELLED', 'UNKNOWN', 'ABORTED')
        def lastStatus = null
        def finalStatus = null

        try {
            spinner.start()

            while (!shouldExit.get() && !Thread.currentThread().isInterrupted()) {
                try {
                    // Fetch workflow status and logs
                    final status = fetchWorkflowStatus(workflowId, queryParams, accessToken, apiEndpoint)
                    final logEntries = fetchWorkflowLogs(workflowId, queryParams, accessToken, apiEndpoint)

                    // Update spinner with status if it changed
                    if (status && status != lastStatus) {
                        def spinnerMode = getSpinnerMode(status, workflowFinished)
                        spinner.updateMessage(formatWorkflowStatus(status), getColorForStatus(status), spinnerMode)
                        lastStatus = status
                    }

                    // Display new log entries
                    def newDisplayedCount = displayLogEntries(logEntries, displayedLogCount, firstLogReceived, shouldExit, spinner)
                    if (newDisplayedCount > displayedLogCount) {
                        if (!firstLogReceived) {
                            firstLogReceived = true
                        }
                        displayedLogCount = newDisplayedCount
                    }

                    // Track workflow completion
                    if (status && finalStatuses.contains(status) && !workflowFinished) {
                        log.debug "Workflow reached final status: ${status}, continuing to poll for ${LOG_GRACE_PERIOD_MS}ms to capture remaining logs"
                        workflowFinished = true
                        workflowFinishedTime = System.currentTimeMillis()
                        finalStatus = status
                    }

                    // Stop after grace period
                    if (workflowFinished && (System.currentTimeMillis() - workflowFinishedTime) > LOG_GRACE_PERIOD_MS) {
                        log.debug "Grace period elapsed, stopping log polling"
                        break
                    }

                    // Wait before next poll
                    sleepInterruptibly(LOG_POLL_INTERVAL_MS, shouldExit)

                }
                catch (IOException e) {
                    log.debug "Log polling got interrupted - ${e.message}"
                    Thread.currentThread().interrupt()
                    break
                }
                catch (Exception e) {
                    if (shouldExit.get())
                        break
                    log.error "Error polling workflow logs: ${e.message}", e
                    sleepInterruptibly(LOG_POLL_INTERVAL_MS, shouldExit)
                }
            }
        } finally {
            if (spinner.isRunning()) {
                spinner.stop()
                // Clear the separator line that was above the spinner
                if (ColorUtil.isAnsiEnabled() && firstLogReceived) {
                    // Move up 2 lines (spinner line + separator line) and clear them
                    print("\u001B[2A\u001B[0J")
                    System.out.flush()
                }
            }

            // Print final status after logs
            if (finalStatus) {
                def colorName = getColorForStatus(finalStatus)
                def coloredIcon = (colorName == 'cyan') ? '⠿' : ColorUtil.colorize('⠿', colorName)
                println "${coloredIcon} ${formatWorkflowStatus(finalStatus)}"
            }
        }
    }

    /**
     * Get the spinner mode for a workflow status
     */
    private String getSpinnerMode(String status, boolean workflowFinished) {
        if (!status) return 'waiting'

        // Use succeeded animation for successful completion
        if (status == 'SUCCEEDED') {
            return 'succeeded'
        }

        // Use failed animation for failed/cancelled/aborted statuses
        if (status == 'FAILED' || status == 'CANCELLED' || status == 'ABORTED' || status == 'UNKNOWN') {
            return 'failed'
        }

        // Use waiting animation for pending/submitted
        if (status == 'PENDING' || status == 'SUBMITTED') {
            return 'waiting'
        }

        // Use running animation for all other active states (RUNNING, etc.)
        return 'running'
    }

    /**
     * Get the color name for a workflow status
     */
    private String getColorForStatus(String status) {
        if (!status) return 'cyan'

        switch (status) {
            case 'PENDING':
            case 'SUBMITTED':
                return 'yellow'
            case 'RUNNING':
                return 'blue'
            case 'FAILED':
            case 'ABORTED':
            case 'CANCELLED':
                return 'red'
            case 'SUCCEEDED':
                return 'green'
            default:
                return 'cyan'
        }
    }

    /**
     * Format workflow status with appropriate color
     */
    private String formatWorkflowStatus(String status) {
        if (!status) return ""

        def colorName = getColorForStatus(status)
        def coloredStatus = (colorName == 'cyan') ? status : ColorUtil.colorize(status, colorName)

        return "Workflow status: ${coloredStatus} "
    }

    /**
     * Fetch workflow status from API
     */
    private String fetchWorkflowStatus(String workflowId, Map queryParams, String accessToken, String apiEndpoint) {
        final workflowResponse = apiGet("/workflow/${workflowId}", queryParams, accessToken, apiEndpoint)
        final workflow = workflowResponse.workflow as Map
        final status = workflow?.status as String
        log.debug "Workflow status: ${status}"
        return status
    }

    /**
     * Fetch workflow logs from API
     */
    private List<String> fetchWorkflowLogs(String workflowId, Map queryParams, String accessToken, String apiEndpoint) {
        final logResponse = apiGet("/workflow/${workflowId}/log", queryParams, accessToken, apiEndpoint)
        final logData = logResponse.log as Map
        return logData?.entries as List<String> ?: []
    }

    /**
     * Display new log entries and return new count
     */
    private int displayLogEntries(List<String> allEntries, int currentCount, boolean firstReceived,
                                   AtomicBoolean shouldExit, SpinnerUtil spinner) {
        if (!allEntries || allEntries.size() <= currentCount) {
            // No new entries - spinner handles the waiting indication
            return currentCount
        }

        // Stop spinner and clear separator line if it was running
        def wasRunning = spinner.isRunning()
        if (wasRunning) {
            spinner.stop()
            // Clear the separator line that was above the spinner
            if (ColorUtil.isAnsiEnabled()) {
                // Move up 2 lines (spinner line + separator line) and clear them
                print("\u001B[2A\u001B[0J")
                System.out.flush()
            }
        }

        // Print separator on first log entry
        if (!firstReceived) {
            println ""
            println ('─' * 70)
            println ""
        }

        // Display new entries
        final newEntries = allEntries.subList(currentCount, allEntries.size())
        for (String entry : newEntries) {
            if (shouldExit.get()) break
            println entry
        }

        // Restart spinner below the logs with separator
        if (wasRunning && !shouldExit.get()) {
            println ""
            println ('─' * 70)
            spinner.start()
        }

        return allEntries.size()
    }

    /**
     * Sleep with ability to be interrupted by exit flag or thread interruption.
     *
     * When Thread.sleep() is interrupted, it clears the interrupt status and throws InterruptedException.
     * We catch it and restore the interrupt status by calling Thread.currentThread().interrupt().
     * This preserved interrupt status is then detected by the outer polling loop at line 506,
     * which checks !Thread.currentThread().isInterrupted() and breaks out gracefully.
     */
    private void sleepInterruptibly(int milliseconds, AtomicBoolean shouldExit) {
        final iterations = milliseconds / LOG_SLEEP_INTERVAL_MS
        try {
            for (int i = 0; i < iterations && !shouldExit.get(); i++) {
                Thread.sleep(LOG_SLEEP_INTERVAL_MS)
            }
        } catch (InterruptedException e) {
            // Restore interrupt status so the outer loop can detect it and exit cleanly
            log.debug "Log sleep interrupted (restoring interrupt status) - ${e.message}"
            Thread.currentThread().interrupt()
        }
    }

    /**
     * Print exit message when log polling is interrupted
     */
    private void printLogPollingExitMessage(String trackingUrl) {
        if (ColorUtil.isAnsiEnabled()) {
            println ""
            println "Exiting log viewer."
            if (trackingUrl) {
                print(ColorUtil.colorize("⚠️ Workflow is still running in Seqera Platform: ", "yellow bold"))
                ColorUtil.printColored(trackingUrl, "cyan")
            }
        } else {
            log.info "Exiting log viewer."
            if (trackingUrl) {
                log.info "Workflow is still running in Seqera Platform: ${trackingUrl}"
            }
        }
    }

    /**
     * Remove shutdown hook if we exit normally
     */
    private void removeShutdownHook(Thread shutdownHook) {
        try {
            Runtime.getRuntime().removeShutdownHook(shutdownHook)
        } catch (IllegalStateException e) {
            // Shutdown hook already executed, ignore
        }
    }

    // ===== Parameter & Config Parsing Methods =====

    /**
     * Build parameters JSON text from CLI params and params file
     */
    private String buildParamsText(Map<String, String> params, String paramsFile) {
        final result = [:]

        // Apply params file
        if (paramsFile) {
            log.debug "Processing params file: ${paramsFile}"
            def path = validateParamsFile(paramsFile)
            def type = path.extension.toLowerCase() ?: null
            if (type == 'json') {
                log.debug "Reading JSON params file: ${paramsFile}"
                readJsonFile(path, result)
            } else if (type == 'yml' || type == 'yaml') {
                log.debug "Reading YAML params file: ${paramsFile}"
                readYamlFile(path, result)
            }
            log.debug "Loaded ${result.size()} parameters from file"
        }

        // Apply CLI params
        if (params) {
            log.debug "Processing ${params.size()} CLI parameters"
            for (Map.Entry<String, String> entry : params) {
                log.trace "Adding parameter: ${entry.key} = ${entry.value}"
                addParam(result, entry.key, entry.value)
            }
            log.debug "Total parameters after CLI merge: ${result.size()}"
        }

        if (result.isEmpty()) {
            log.debug "No parameters to include in launch request"
            return null
        }

        // Convert to JSON
        def jsonBuilder = new JsonBuilder(result)
        def jsonText = jsonBuilder.toString()
        log.debug "Converted ${result.size()} parameters to JSON format"
        return jsonText
    }

    private String buildConfigText(List<String> configFiles) {
        if (!configFiles || configFiles.isEmpty()) {
            log.debug "No configuration files specified"
            return null
        }

        log.debug "Processing ${configFiles.size()} configuration files"
        def configText = new StringBuilder()
        for (String configFile : configFiles) {
            log.debug "Reading config file: ${configFile}"
            def path = FileHelper.asPath(configFile)
            if (!path.exists()) {
                log.debug "Config file not found: ${configFile}"
                throw new AbortOperationException("Config file not found: ${configFile}")
            }
            def content = path.text
            log.debug "Read ${content.length()} characters from config file: ${configFile}"
            configText.append(content).append('\n')
        }

        log.debug "Combined configuration text: ${configText.length()} characters"
        return configText.toString()
    }

    private Path validateParamsFile(String file) {
        def result = FileHelper.asPath(file)
        def ext = result.getExtension()
        if (!VALID_PARAMS_FILE.contains(ext)) {
            throw new AbortOperationException("Not a valid params file extension: $file -- It must be one of the following: ${VALID_PARAMS_FILE.join(',')}")
        }
        return result
    }

    private void readJsonFile(Path file, Map result) {
        try {
            def json = (Map<String, Object>) new JsonSlurper().parseText(file.text)
            json.forEach((name, value) -> {
                addParam0(result, name, value)
            })
        } catch (NoSuchFileException | FileNotFoundException e) {
            throw new AbortOperationException("Specified params file does not exist: ${file.toUriString()}")
        } catch (Exception e) {
            throw new AbortOperationException("Cannot parse params file: ${file.toUriString()} - Cause: ${e.message}", e)
        }
    }

    private void readYamlFile(Path file, Map result) {
        try {
            def yaml = (Map<String, Object>) new Yaml().load(file.text)
            yaml.forEach((name, value) -> {
                addParam0(result, name, value)
            })
        } catch (NoSuchFileException | FileNotFoundException e) {
            throw new AbortOperationException("Specified params file does not exist: ${file.toUriString()}")
        } catch (Exception e) {
            throw new AbortOperationException("Cannot parse params file: ${file.toUriString()}", e)
        }
    }

    private static void addParam(Map params, String key, String value, List path = [], String fullKey = null) {
        if (!fullKey) {
            fullKey = key
        }
        final m = DOT_NOT_ESCAPED.matcher(key)
        if (m.find()) {
            final p = m.start()
            final root = key.substring(0, p)
            if (!root) throw new AbortOperationException("Invalid parameter name: $fullKey")
            path.add(root)
            def nested = params.get(root)
            if (nested == null) {
                nested = new LinkedHashMap<>()
                params.put(root, nested)
            } else if (!(nested instanceof Map)) {
                log.warn "Command line parameter --${path.join('.')} is overwritten by --${fullKey}"
                nested = new LinkedHashMap<>()
                params.put(root, nested)
            }
            addParam((Map) nested, key.substring(p + 1), value, path, fullKey)
        } else {
            addParam0(params, key.replaceAll(DOT_ESCAPED, '.'), parseParamValue(value))
        }
    }

    private static void addParam0(Map params, String key, Object value) {
        if (key.contains('-')) {
            key = kebabToCamelCase(key)
        }
        params.put(key, value)
    }

    private static String kebabToCamelCase(String str) {
        final result = new StringBuilder()
        str.split('-').eachWithIndex { String entry, int i ->
            result << (i > 0 ? entry.capitalize() : entry)
        }
        return result.toString()
    }

    private static Object parseParamValue(String str) {
        if (str == null) return null

        if (str.toLowerCase() == 'true') return Boolean.TRUE
        if (str.toLowerCase() == 'false') return Boolean.FALSE

        if (str ==~ /-?\d+(\.\d+)?/ && str.isInteger()) return str.toInteger()
        if (str ==~ /-?\d+(\.\d+)?/ && str.isLong()) return str.toLong()
        if (str ==~ /-?\d+(\.\d+)?/ && str.isDouble()) return str.toDouble()

        return str
    }

    // ===== Utility Methods =====

    /**
     * Check if a path is a local file path (vs remote repository)
     */
    private static boolean isLocalPath(String path) {
        if (path.startsWith('/') || path.startsWith('./') || path.startsWith('../')) {
            return true
        }
        if (path.contains('\\') || path ==~ /^[A-Za-z]:.*/) {
            return true
        }
        return false
    }

    /**
     * Resolve pipeline name to full repository URL using AssetManager
     */
    protected String resolvePipelineUrl(String pipelineName) {
        try {
            log.debug "Resolving pipeline name using AssetManager: ${pipelineName}"
            def assetManager = new AssetManager(pipelineName)
            def repositoryUrl = assetManager.getRepositoryUrl()
            if (repositoryUrl) {
                log.debug "AssetManager resolved URL: ${repositoryUrl}"
                return repositoryUrl
            } else {
                log.debug "AssetManager could not resolve URL, using original name: ${pipelineName}"
                return pipelineName
            }
        } catch (Exception e) {
            log.debug "Failed to resolve pipeline URL with AssetManager: ${e.message}, using original name: ${pipelineName}"
            return pipelineName
        }
    }

    // ===== API Helper Methods =====

    protected Map apiGet(String path, Map queryParams = [:], String accessToken, String apiEndpoint) {
        final url = buildUrl(apiEndpoint, path, queryParams)
        log.debug "Platform API - GET ${url}"
        final client = createHttpClient(accessToken)
        final request = HttpRequest.newBuilder()
            .uri(URI.create(url))
            .GET()
            .build()

        final response = client.send(request, HttpResponse.BodyHandlers.ofString())

        if (response.statusCode() != 200) {
            final error = response.body() ?: "HTTP ${response.statusCode()}"
            throw new RuntimeException("API request failed: ${error}")
        }

        return new JsonSlurper().parseText(response.body()) as Map
    }

    protected Map apiPost(String path, Map body, Map queryParams = [:], String accessToken, String apiEndpoint) {
        final url = buildUrl(apiEndpoint, path, queryParams)
        log.debug "Platform API - POST ${url}"
        final requestBody = new JsonBuilder(body).toString()
        final client = createHttpClient(accessToken)
        final request = HttpRequest.newBuilder()
            .uri(URI.create(url))
            .header('Content-Type', 'application/json')
            .POST(HttpRequest.BodyPublishers.ofString(requestBody))
            .build()

        final response = client.send(request, HttpResponse.BodyHandlers.ofString())

        if (response.statusCode() != 200) {
            if (response.statusCode() == 403) {
                throw new AbortOperationException("ERROR: Unable to launch workflow.\nCheck your credentials with 'nextflow auth status' and your user role in the workspace (required: 'maintain' or higher).")
            }
            final error = response.body() ?: "HTTP ${response.statusCode()}"
            throw new RuntimeException("Failed to launch workflow: ${error}")
        }

        return new JsonSlurper().parseText(response.body()) as Map
    }

    private String buildUrl(String endpoint, String path, Map queryParams) {
        def url = new StringBuilder(endpoint)
        if (!path.startsWith('/')) {
            url.append('/')
        }
        url.append(path)

        if (queryParams && !queryParams.isEmpty()) {
            url.append('?')
            url.append(queryParams.collect { k, v -> "${URLEncoder.encode(k.toString(), 'UTF-8')}=${URLEncoder.encode(v.toString(), 'UTF-8')}" }.join('&'))
        }

        return url.toString()
    }

    // ===== Workspace & User Helper Methods =====

    protected Long resolveWorkspaceId(Map config, String workspaceName, String accessToken, String apiEndpoint) {
        // First check config for workspace ID
        final configWorkspaceId = config['tower.workspaceId']
        if (configWorkspaceId) {
            return configWorkspaceId as Long
        }

        // If workspace name provided, look it up
        if (workspaceName) {
            final userId = getUserInfo(accessToken, apiEndpoint).id as String
            final workspaces = listUserWorkspaces(accessToken, apiEndpoint, userId)

            final matchingWorkspace = workspaces.find { workspace ->
                final ws = workspace as Map
                final wsName = ws.workspaceName as String
                wsName == workspaceName
            }

            if (!matchingWorkspace) {
                throw new AbortOperationException("Workspace '${workspaceName}' not found")
            }

            return (matchingWorkspace as Map).workspaceId as Long
        }

        return null
    }

    protected Map findComputeEnv(String computeEnvName, Long workspaceId, String accessToken, String apiEndpoint) {
        final computeEnvs = listComputeEnvironments( accessToken, apiEndpoint, workspaceId ? workspaceId.toString() : null)

        log.debug "Looking for ${computeEnvName ? "compute environment with name: ${computeEnvName}" : "primary compute environment"} ${workspaceId ? "in workspace ID ${workspaceId}" : "in personal workspace"}"

        for (item in computeEnvs) {
            final computeEnv = item as Map
            if ((computeEnvName && computeEnv.name == computeEnvName) || (!computeEnvName && computeEnv.primary == true)) {
                log.debug "Found ${computeEnvName ? "matching" : "primary"} compute environment '${computeEnv.name}' (ID: ${computeEnv.id})"
                return computeEnv
            }
        }

        return null
    }

}
