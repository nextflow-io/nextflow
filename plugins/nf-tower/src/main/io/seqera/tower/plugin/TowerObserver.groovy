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

package io.seqera.tower.plugin

import groovy.transform.CompileStatic
import groovy.transform.ToString
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.container.resolver.ContainerMeta
import nextflow.exception.AbortOperationException
import nextflow.processor.TaskHandler
import nextflow.processor.TaskId
import nextflow.processor.TaskProcessor
import nextflow.script.PlatformMetadata
import nextflow.trace.ResourcesAggregator
import nextflow.trace.TraceObserverV2
import nextflow.trace.TraceRecord
import nextflow.trace.event.FilePublishEvent
import nextflow.trace.event.TaskEvent
import nextflow.util.Duration
import nextflow.util.LoggerHelper
import nextflow.util.ProcessHelper
import nextflow.util.Threads

import java.time.Instant
import java.time.OffsetDateTime
import java.time.ZoneId
import java.util.concurrent.ConcurrentHashMap
import java.util.concurrent.LinkedBlockingQueue
import java.util.concurrent.TimeUnit

/**
 * Send out messages via HTTP to a configured URL on different workflow
 * execution events.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class TowerObserver implements TraceObserverV2 {

    static private final int TASKS_PER_REQUEST = 100

    static private final Duration REQUEST_INTERVAL = Duration.of('1 sec')

    static private final Duration ALIVE_INTERVAL = Duration.of('1 min')

    @ToString(includeNames = true)
    static class ProcessEvent {
        TraceRecord trace
        boolean completed
    }

    private Session session

    /**
     * Workflow identifier, will be taken from the Session() object later
     */
    private String runName

    /**
     * Store the sessions unique ID for downstream reference purposes
     */
    private String runId

    private String workflowId

    private String watchUrl

    private ResourcesAggregator aggregator

    protected Map<String,String> env = System.getenv()

    private LinkedBlockingQueue<ProcessEvent> events = new LinkedBlockingQueue()

    private Thread sender

    private Duration requestInterval = REQUEST_INTERVAL

    private Duration aliveInterval = ALIVE_INTERVAL

    private LinkedHashSet<String> processNames = new LinkedHashSet<>(20)

    private boolean towerLaunch

    private String workspaceId

    private TowerReports reports

    private TowerClient client

    private Map<String,Boolean> allContainers = new ConcurrentHashMap<>()

    TowerObserver(Session session, TowerClient client, String workspaceId, Map env) {
        this.session = session
        this.workspaceId = workspaceId
        this.reports = new TowerReports(session)
        this.client = client
        if( env )
            this.env = env
    }


    @Override
    boolean enableMetrics() { true }

    String getWorkflowId() { workflowId }

    boolean getTowerLaunch() { towerLaunch }

    String getRunName() { runName }

    String getRunId() { runId }

    void setAliveInterval(Duration d) {
        this.aliveInterval = d
    }

    void setRequestInterval(Duration d) {
        this.requestInterval = d
    }

    String getWorkspaceId() { workspaceId }


    /**
     * On workflow start, submit a message with some basic
     * information, like Id, activity and an ISO 8601 formatted
     * timestamp.
     * @param session The current Nextflow session object
     */
    @Override
    void onFlowCreate(Session session) {
        log.debug "Creating Seqera Platform observer -- endpoint=$client.endpoint; requestInterval=$requestInterval; aliveInterval=$aliveInterval;"

        this.session = session
        this.aggregator = new ResourcesAggregator()
        this.runName = session.getRunName()
        this.runId = session.getUniqueId()

        // send hello to verify auth
        final ret = client.traceCreate(makeCreateReq(session), workspaceId)
        this.workflowId = ret.workflowId
        if( !workflowId )
            throw new AbortOperationException("Invalid Seqera Platform API response - Missing workflow Id")
        log.debug "Platform workflow id: $workflowId; workflow url: ${ret.watchUrl}"
        session.workflowMetadata.platform.workflowId = workflowId
        // note: `watchUrl` in the create response requires Platform 26.01 or later
        this.watchUrl = ret.watchUrl as String
        session.workflowMetadata.platform.workflowUrl = watchUrl
        if( ret.message )
            log.warn(ret.message.toString())
        // populate platform metadata from the create response
        if( ret.metadata )
            applyPlatformMetadata(ret.metadata as Map)

        // Prepare to collect report paths if tower configuration has a 'reports' section
        reports.flowCreate(workflowId)
    }


    /**
     * Apply platform metadata received inline from the trace create response.
     * This avoids extra API calls to fetch user, workspace, and launch details.
     */
    protected void applyPlatformMetadata(Map metadata) {
        try {
            final platform = session.workflowMetadata.platform
            // user info
            if( metadata.userId )
                platform.user = new PlatformMetadata.User(
                    id: metadata.userId as String,
                    userName: metadata.userName as String,
                    organization: metadata.userOrganization as String
                )
            // workspace info
            if( metadata.workspaceId )
                platform.workspace = new PlatformMetadata.Workspace(
                    workspaceId: metadata.workspaceId as String,
                    workspaceName: metadata.workspaceName as String,
                    workspaceFullName: metadata.workspaceFullName as String,
                    orgName: metadata.orgName as String
                )
            // launch details (only present for Platform-submitted runs)
            if( metadata.computeEnvId )
                platform.computeEnv = new PlatformMetadata.ComputeEnv(
                    id: metadata.computeEnvId as String,
                    name: metadata.computeEnvName as String,
                    platform: metadata.computeEnvPlatform as String
                )
            if( metadata.pipelineName )
                platform.pipeline = new PlatformMetadata.Pipeline(
                    id: metadata.pipelineId as String,
                    name: metadata.pipelineName as String,
                    revision: metadata.revision as String,
                    commitId: metadata.commitId as String
                )
            if( metadata.labels )
                platform.labels = metadata.labels as List<String>
        }
        catch( Exception e ) {
            log.debug("Failed to apply platform metadata from create response", e)
        }
    }

    protected Map makeCreateReq(Session session) {
        def result = new HashMap(5)
        result.sessionId = session.uniqueId.toString()
        result.runName = session.runName
        result.projectName = session.workflowMetadata.projectName
        result.repository = session.workflowMetadata.repository
        result.workflowId = env.get('TOWER_WORKFLOW_ID')
        result.instant = Instant.now().toEpochMilli()
        this.towerLaunch = result.workflowId != null
        return result
    }

    @Override
    void onProcessCreate(TaskProcessor process) {
        log.trace "Creating process ${process.name}"
        if( !processNames.add(process.name) )
            throw new IllegalStateException("Process name `${process.name}` already used")
    }

    @Override
    void onFlowBegin() {
        // configure error retry

        final payload = client.traceBegin(makeBeginReq(session), workspaceId, workflowId)
        this.watchUrl ?= payload.watchUrl
        session.workflowMetadata.platform.workflowUrl ?= watchUrl
        this.sender = Threads.start('Tower-thread', this.&sendTasks0)
        final msg = "Monitor the execution with Seqera Platform using this URL: ${watchUrl}"
        log.info(LoggerHelper.STICKY, msg)
    }

    String getAccessToken() {
        if( !accessToken )
            throw new AbortOperationException("Missing Seqera Platform access token -- Make sure there's a variable TOWER_ACCESS_TOKEN in your environment")
        return accessToken
    }

    /**
     * Send an HTTP message when the workflow is completed.
     */
    @Override
    void onFlowComplete() {
        // publish runtime reports
        reports.publishRuntimeReports()
        // submit the completion record
        if( sender ) {
            events << new ProcessEvent(completed: true)
            // wait the submission of pending events
            sender.join()
        }
        // wait and flush reports content
        reports.flowComplete()
        // notify the workflow completion
        // note: only send complete if onFlowBegin was invoked (sender is set there)
        if( workflowId && sender ) {
            client.traceComplete(makeCompleteReq(session), workspaceId, workflowId)
        }
    }

    @Override
    void onTaskPending(TaskEvent event) {
        events << new ProcessEvent(trace: event.trace)
    }

    /**
     * Send an HTTP message when a process has been submitted
     *
     * @param handler A {@link TaskHandler} object representing the task submitted
     * @param trace A {@link TraceRecord} object holding the task metadata and runtime info
     */
    @Override
    void onTaskSubmit(TaskEvent event) {
        events << new ProcessEvent(trace: event.trace)
    }

    /**
     * Send an HTTP message, when a process has started
     *
     * @param handler A {@link TaskHandler} object representing the task started
     * @param trace A {@link TraceRecord} object holding the task metadata and runtime info
     */
    @Override
    void onTaskStart(TaskEvent event) {
        events << new ProcessEvent(trace: event.trace)
    }

    /**
     * Send an HTTP message, when a process completed
     *
     * @param handler A {@link TaskHandler} object representing the task completed
     * @param trace A {@link TraceRecord} object holding the task metadata and runtime info
     */
    @Override
    void onTaskComplete(TaskEvent event) {
        events << new ProcessEvent(trace: event.trace)

        synchronized (this) {
            aggregator.aggregate(event.trace)
        }
    }

    @Override
    void onTaskCached(TaskEvent event) {
        // event was triggered by a stored task, ignore it
        if( !event.trace )
            return

        // add the cached task event
        events << new ProcessEvent(trace: event.trace)

        // remove the record from the current records
        synchronized (this) {
            aggregator.aggregate(event.trace)
        }
    }

    /**
     * Send an HTTP message, when a workflow has failed
     *
     * @param handler A {@link TaskHandler} object representing the task that caused the workflow execution to fail (it may be null)
     * @param trace A {@link TraceRecord} object holding the task metadata and runtime info (it may be null)
     */
    @Override
    void onFlowError(TaskEvent event) {
        events << new ProcessEvent(trace: event.trace)
    }

    /**
     * Update reports file when a file is published
     *
     * @param destination File path at `publishDir` of the published file.
     */
    @Override
    void onFilePublish(FilePublishEvent event) {
        reports.filePublish(event.target)
    }

    protected boolean isCliLogsEnabled() {
        return env.get('TOWER_ALLOW_NEXTFLOW_LOGS') == 'true'
    }

    protected String getOperationId() {
        if( !isCliLogsEnabled() )
            return null
        try {
            if( env.get('AWS_BATCH_JOB_ID') )
                return  "aws-batch::${env.get('AWS_BATCH_JOB_ID')}"
            else
                return "local-platform::${ProcessHelper.selfPid()}"
        }
        catch (Exception e) {
            log.warn "Unable to retrieve native environment operation id", e
            return null
        }
    }

    protected String getLogFile() {
        return isCliLogsEnabled() ? env.get('NXF_LOG_FILE') : null
    }

    protected String getOutFile() {
        return isCliLogsEnabled() ? env.get('NXF_OUT_FILE') : null
    }

    protected Map makeBeginReq(Session session) {
        def workflow = session.getWorkflowMetadata().toMap()
        workflow.params = session.getParams()
        workflow.id = getWorkflowId()
        workflow.remove('stats')

        // render as a string
        workflow.container = mapToString(workflow.container)
        workflow.configText = session.resolvedConfig
        // extra metadata
        workflow.operationId = getOperationId()
        workflow.logFile = getLogFile()
        workflow.outFile = getOutFile()

        def result = new LinkedHashMap(5)
        result.workflow = workflow
        result.processNames = new ArrayList(processNames)
        result.towerLaunch = towerLaunch
        result.instant = Instant.now().toEpochMilli()
        return result
    }

    protected Map makeCompleteReq(Session session) {
        def workflow = session.getWorkflowMetadata().toMap()
        //Remove retrieved platform info
        if( workflow.platform )
            workflow.remove('platform')

        workflow.params = session.getParams()
        workflow.id = getWorkflowId()
        // render as a string
        workflow.container = mapToString(workflow.container)
        workflow.configText = session.resolvedConfig
        // extra metadata
        workflow.operationId = getOperationId()
        workflow.logFile = getLogFile()
        workflow.outFile = getOutFile()

        def result = new LinkedHashMap(5)
        result.workflow = workflow
        result.metrics = getMetricsList()
        result.progress = getWorkflowProgress(false)
        result.instant = Instant.now().toEpochMilli()
        return result
    }

    protected Map makeHeartbeatReq() {
        def result = new HashMap(1)
        result.progress = getWorkflowProgress(true)
        result.instant = Instant.now().toEpochMilli()
        return result
    }

    protected String mapToString(def obj) {
        if( obj == null )
            return null
        if( obj instanceof CharSequence )
            return obj.toString()
        if( obj instanceof Map ) {
            // turn this off for multiple containers because the string representation is broken
            return null
        }
        throw new IllegalArgumentException("Illegal container attribute type: ${obj.getClass().getName()} = ${obj}" )
    }

    protected Map makeTaskMap0(TraceRecord trace) {
        Map<String,?> record = new LinkedHashMap<>(trace.store.size())
        for( Map.Entry<String,Object> entry : trace.store.entrySet() ) {
            def name = entry.key
            // remove '%' char from field prefix
            if( name.startsWith('%') )
                name = 'p' + name.substring(1)
            // normalise to camelCase
            name = underscoreToCamelCase(name)
            // put the value
            record.put(name, fixTaskField(name,entry.value))
        }

        // prevent invalid tag data
        if( record.tag!=null && !(record.tag instanceof CharSequence)) {
            final msg = "Invalid tag value for process: ${record.process} -- A string is expected instead of type: ${record.tag.getClass().getName()}; offending value=${record.tag}"
            log.warn1(msg, cacheKey: record.process)
            record.tag = null
        }

        // add transient fields
        record.executor = trace.getExecutorName()
        record.cloudZone = trace.getMachineInfo()?.zone
        record.machineType = trace.getMachineInfo()?.type
        record.priceModel = trace.getMachineInfo()?.priceModel?.toString()
        record.numSpotInterruptions = trace.getNumSpotInterruptions()
        record.logStreamId = trace.getLogStreamId()

        return record
    }


    static protected Object fixTaskField(String name, value) {
        if( TraceRecord.FIELDS[name] == 'date' )
            return value ? OffsetDateTime.ofInstant(Instant.ofEpochMilli(value as long), ZoneId.systemDefault()) : null
        else
            return value
    }

    protected Map makeTasksReq(Collection<TraceRecord> tasks) {

        def payload = new ArrayList(tasks.size())
        for( TraceRecord rec : tasks ) {
            payload << makeTaskMap0(rec)
        }

        final result = new LinkedHashMap(5)
        result.put('tasks', payload)
        result.put('progress', getWorkflowProgress(true))
        result.put('containers', getNewContainers(tasks))
        result.instant = Instant.now().toEpochMilli()
        return result
    }

    protected List<ContainerMeta> getNewContainers(Collection<TraceRecord> tasks) {
        final result = new ArrayList<ContainerMeta>()
        for( TraceRecord it : tasks ) {
            final meta = it.getContainerMeta()
            if( meta && !allContainers.get(meta.targetImage) ) {
                allContainers.put(meta.targetImage, Boolean.TRUE)
                result.add(meta)
            }
        }
        return result
    }

    protected List getMetricsList() {
        return aggregator.computeSummaryList()
    }

    protected WorkflowProgress getWorkflowProgress(boolean quick) {
        def stats = quick ? session.getStatsObserver().getQuickStats() : session.getStatsObserver().getStats()
        new WorkflowProgress(stats)
    }

    protected String underscoreToCamelCase(String str) {
        if( !str.contains('_') )
            return str

        final words = str.tokenize('_')
        def result = words[0]
        for( int i=1; i<words.size(); i++ )
            result+=words[i].capitalize()

        return result
    }


    protected Map<String,Integer> loadSchema() {
        final props = new Properties()
        props.load(this.getClass().getResourceAsStream('/tower-schema.properties'))
        final result = new HashMap<String,Integer>(props.size())
        for( String key : props.keySet() ) {
            final value = props.getProperty(key)
            result.put( key, value ? value as Integer : null )
        }
        return result
    }

    protected void sendTasks0(dummy) {
        final tasks = new HashMap<TaskId, TraceRecord>(TASKS_PER_REQUEST)
        boolean complete = false
        long previous = System.currentTimeMillis()
        final long period = requestInterval.millis
        final long delay = period / 10 as long

        while( !complete ) {
            final ProcessEvent ev = events.poll(delay, TimeUnit.MILLISECONDS)
            // reconcile task events ie. send out only the last event
            if( ev ) {
                log.trace "Tower event=$ev"
                if( ev.trace )
                    tasks[ev.trace.taskId] = ev.trace
                if( ev.completed )
                    complete = true
            }

            // check if there's something to send
            final now = System.currentTimeMillis()
            final delta = now -previous

            if( !tasks ) {
                if( delta > aliveInterval.millis ) {
                    final req = makeHeartbeatReq()
                    client.traceHeartbeat(req, workspaceId, workflowId)
                    previous = now
                }
                continue
            }

            if( delta > period || tasks.size() >= TASKS_PER_REQUEST || complete ) {
                // send
                final req = makeTasksReq(tasks.values())
                client.traceProgress(req, workspaceId, workflowId)

                // clean up for next iteration
                previous = now
                tasks.clear()
            }
        }
    }

}
