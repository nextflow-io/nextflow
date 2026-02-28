/*
 * Copyright 2013-2026, Seqera Labs
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

package nextflow.trace

import java.util.concurrent.ConcurrentHashMap

import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.processor.TaskRun
import nextflow.trace.event.TaskEvent

/**
 * AI agent-friendly log observer that outputs minimal, structured information
 * optimized for AI context windows.
 *
 * Activated via environment variables: NXF_AGENT=1, AGENT=1, or CLAUDECODE=1
 *
 * Output format:
 * - [PIPELINE] name version | profile=X
 * - [WORKDIR] /path/to/work
 * - [WARN] warning message (deduplicated)
 * - [ERROR] PROCESS_NAME (sample_id) with exit/cmd/stderr/workdir
 * - [SUCCESS|FAILED] completed=N failed=N cached=N
 *
 * @author Edmund Miller <edmund.miller@utdallas.edu>
 */
@Slf4j
class AgentLogObserver implements TraceObserverV2 {

    private Session session
    private WorkflowStatsObserver statsObserver
    private final Set<String> seenWarnings = ConcurrentHashMap.newKeySet()
    private volatile boolean started = false
    private volatile boolean completed = false

    /**
     * Set the workflow stats observer for retrieving task statistics
     */
    void setStatsObserver(WorkflowStatsObserver observer) {
        this.statsObserver = observer
    }

    /**
     * Print a line to standard output (agent format)
     */
    protected void println(String line) {
        System.err.println(line)
    }

    // -- TraceObserverV2 lifecycle methods --

    @Override
    void onFlowCreate(Session session) {
        this.session = session
    }

    @Override
    void onFlowBegin() {
        if( started )
            return
        started = true

        // Print pipeline info
        def manifest = session.manifest
        def pipelineName = manifest?.name ?: session.scriptName ?: 'unknown'
        def version = manifest?.version ?: ''
        def profile = session.profile ?: 'standard'

        def info = "[PIPELINE] ${pipelineName}"
        if( version )
            info += " ${version}"
        info += " | profile=${profile}"
        println(info)

        // Print work directory
        def workDir = session.workDir?.toUriString() ?: session.workDir?.toString()
        if( workDir )
            println("[WORKDIR] ${workDir}")
    }

    @Override
    void onFlowComplete() {
        if( completed )
            return
        completed = true
        printSummary()
    }

    @Override
    void onFlowError(TaskEvent event) {
        // Error is already reported by onTaskComplete for failed tasks
    }

    @Override
    void onTaskSubmit(TaskEvent event) {
        // Not reported in agent mode
    }

    @Override
    void onTaskComplete(TaskEvent event) {
        def handler = event.handler
        def task = handler?.task
        if( task?.isFailed() ) {
            printTaskError(task)
        }
    }

    @Override
    void onTaskCached(TaskEvent event) {
        // Not reported in agent mode
    }

    /**
     * Append a warning message (deduplicated)
     */
    void appendWarning(String message) {
        if( message == null )
            return
        // Normalize and deduplicate
        def normalized = message.trim().replaceAll(/\s+/, ' ')
        if( seenWarnings.add(normalized) ) {
            println("[WARN] ${normalized}")
        }
    }

    /**
     * Append an error message
     */
    void appendError(String message) {
        if( message )
            println("[ERROR] ${message}")
    }

    /**
     * Append info (suppressed in agent mode except for critical messages)
     */
    void appendInfo(String message) {
        // Suppress most info messages in agent mode
        // Only pass through if it contains critical keywords
        if( message && (message.contains('ERROR') || message.contains('WARN')) ) {
            println(message)
        }
    }

    /**
     * Print task error with full diagnostic context
     */
    protected void printTaskError(TaskRun task) {
        def name = task.getName()
        println("[ERROR] ${name}")

        // Exit status
        def exitStatus = task.getExitStatus()
        if( exitStatus != null && exitStatus != Integer.MAX_VALUE ) {
            println("exit: ${exitStatus}")
        }

        // Command/script (first line or truncated)
        def script = task.getScript()?.toString()?.trim()
        if( script ) {
            // Truncate long commands
            def cmd = script.length() > 200 ? script.substring(0, 200) + '...' : script
            cmd = cmd.replaceAll(/\n/, ' ').replaceAll(/\s+/, ' ')
            println("cmd: ${cmd}")
        }

        // Stderr
        def stderr = task.getStderr()
        if( stderr ) {
            def lines = stderr.readLines()
            if( lines.size() > 10 ) {
                lines = lines[-10..-1]
            }
            println("stderr: ${lines.join(' | ')}")
        }

        // Stdout (only if relevant)
        def stdout = task.getStdout()
        if( stdout && !stderr ) {
            def lines = stdout.readLines()
            if( lines.size() > 5 ) {
                lines = lines[-5..-1]
            }
            println("stdout: ${lines.join(' | ')}")
        }

        // Work directory
        def workDir = task.getWorkDir()
        if( workDir ) {
            println("workdir: ${workDir.toUriString()}")
        }
    }

    /**
     * Print final summary line
     */
    protected void printSummary() {
        def stats = statsObserver?.getStats()
        def succeeded = stats?.succeededCount ?: 0
        def failed = stats?.failedCount ?: 0
        def cached = stats?.cachedCount ?: 0
        def completed = succeeded + failed

        def status = failed > 0 ? 'FAILED' : 'SUCCESS'
        println("\n[${status}] completed=${completed} failed=${failed} cached=${cached}")
    }

    /**
     * Force termination - called on abort
     */
    void forceTermination() {
        if( !completed ) {
            completed = true
            printSummary()
        }
    }
}
