/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

package nextflow.scheduler
import java.util.concurrent.ConcurrentLinkedQueue
import java.util.concurrent.Executors
import java.util.concurrent.ScheduledExecutorService
import java.util.concurrent.TimeUnit

import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.cloud.CloudConfig
import nextflow.cloud.CloudDriver
import nextflow.cloud.CloudDriverFactory
import nextflow.cloud.types.CloudInstanceStatus
import nextflow.executor.IgBaseTask
import nextflow.processor.TaskId
import nextflow.scheduler.Protocol.NodeData
import nextflow.scheduler.Protocol.TaskHolder
import nextflow.util.Duration
import org.apache.ignite.Ignite
/**
 * Implements an auto-scaling policy which adds new instances
 * under scheduler request pressure and remove when those
 * instances are idle
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class Autoscaler implements Closeable {

    /**
     * The cloud driver name e.g. `aws`
     */
    private String driverName

    /**
     * The reference to the cloud driver
     */
    private CloudDriver driver

    /**
     * The reference to the underlying Ignite cluster
     */
    private Ignite ignite

    /**
     * Map each node UUID to the associated host name
     */
    private Map<UUID, NodeData> workerNodes

    /**
     * Map of scheduled tasks not yet completed
     */
    private Map<TaskId,TaskHolder> scheduledTasks

    private ScheduledExecutorService watchdog

    /**
     * The auto-scaler configuration object
     */
    private CloudConfig.Autoscale scalerConfig

    /**
     * Whenever the auto-scaling is enabled
     */
    private volatile boolean enabled

    /**
     * The list of remaining instances requested by the auto-scaler
     * pending to join the cluster
     */
    private volatile Queue<String> pendingInstanceIds

    private Map<List,Byte> canRunTaskCache = new HashMap<>()

    private long waitingTimestamp

    /**
     * Creates the auto-scaler object
     *
     * @param ignite Reference to the underlying {@link Ignite} cluster
     * @param config A {@link CloudConfig} object holding the cloud configuration
     */
    Autoscaler(Ignite ignite, CloudConfig config) {
        this.ignite = ignite
        this.scalerConfig = config.getAutoscale()
        this.driverName = config.getDriverName()
        this.enabled = scalerConfig.enabled
        log.debug "### Auto-scaling enabled: $enabled"
    }

    /**
     * ONLY FOR TESTING PURPOSE
     */
    protected Autoscaler() { }

    /**
     * Initialize the auto-scaling policy
     *
     * @param nodes The map holding the current cluster topology
     * @param tasks The map of scheduled tasks
     */
    @PackageScope
    void init(Map<UUID,NodeData> nodes, Map<TaskId,TaskHolder> tasks) {

        this.workerNodes = nodes
        this.scheduledTasks = tasks

        // -- init the cloud driver
        if( enabled ) {
            def driverName = driverName ?: CloudDriverFactory.getDefaultDriverName()
            if( !driverName )
                throw new IllegalStateException("No cloud driver name has been specified")

            driver = CloudDriverFactory.getDriver(driverName)
            if( !driver ) throw new IllegalStateException("Cannot load cloud driver: `$driverName`")
            driver.validate(scalerConfig)
        }

        // --
        this.watchdog = Executors.newScheduledThreadPool(1)
        this.watchdog.scheduleWithFixedDelay(this.&clusterWatchdog as Runnable, 1, 1, TimeUnit.MINUTES)
    }

    /**
     * Method invoked when a new node join the cluster. If it's an
     * instance requested by the auto-scaler, remove its ID from the list on
     * pending instances to be launched. Once no more instances and pending
     * it re-enable the auto-scaler
     *
     * @param data A {@link NodeData} object representing the new node
     */
    @PackageScope
    void onNodeStart(NodeData data) {
        if( !pendingInstanceIds )
            // nothing to do
            return

        def found = pendingInstanceIds.remove(data.instanceId)
        if( found ) {
            log.debug "### Autoscale node joined the cluster [${data.hostName}] instance-id=$data.instanceId -- Still missing ${pendingInstanceIds.size()} instances"
            if( pendingInstanceIds.isEmpty() ) {
                enabled = true
            }
        }
    }

    /**
     * Implements the auto-scaler main tasks. This method is invoked periodically every minute
     */
    @PackageScope
    void clusterWatchdog() {
        try {
            checkClusterSize()
            checkStarvingTasks()
            checkIdleNodes()
            checkPendingInstances()
        }
        catch (InterruptedException e) {
            log.debug "Shutdown in progress [InterruptedException]"
        }
        catch( Throwable e ) {
            log.debug "### Oops.. Something went wrong", e
        }
    }

    /**
     * @return the current number of nodes that made-up the cloud cluster
     */
    @PackageScope
    int getClusterSize() {
        return workerNodes.size()
    }

    /**
     * Check that the cluster size is within its min - max size boundaries
     */
    @PackageScope
    void checkClusterSize() {
        if( !enabled ) return

        def currentSize = getClusterSize()
        if( currentSize < scalerConfig.minInstances ) {
            // -- add instances as needed to reach the `minInstances` count
            def missingInstances = scalerConfig.minInstances - currentSize
            log.debug "### Adding $missingInstances instance(s) to cloud cluster -- current-size: $currentSize; min-size: ${scalerConfig.minInstances}; "
            requestNewInstances(missingInstances)
        }
    }

    /**
     * Launch new cloud instances if tasks are not executed after a specified amount of time
     */
    @PackageScope
    void checkStarvingTasks() {

        // find all waiting tasks i.e. not marked as `started`
        final timeout = scalerConfig.starvingTimeout
        final type = scalerConfig.instanceType
        Collection<TaskHolder> waiting = scheduledTasks.values().findAll{ !it.started }

        // no task waiting -- reset the timestamp
        if( !waiting ) {
            log.trace "### No tasks waiting for execution"
            waitingTimestamp = 0
            return
        }

        // some tasks are waiting -- set the timestamp
        if( !waitingTimestamp ) {
            waitingTimestamp = System.currentTimeMillis()
        }

        def delta = System.currentTimeMillis() - waitingTimestamp
        log.trace "### Tasks waiting for execution: count=${waiting.size()}; delta=${Duration.of(delta)}; timeout=${timeout}"
        if( delta < timeout.millis ) {
            return
        }

        // make sure these tasks can run on the instance type chosen
        // and count the number of required cpus
        List<TaskId> tasks = []
        int missingCpus = 0
        def itr = waiting.iterator()
        while( itr.hasNext() ) {
            def it = itr.next()
            if( it.isWaitingMoreThan(timeout) )
                canRunOnExistingNodes(it.task)

            def result = enabled && canRunOnNewInstance(it.task, type)
            if( result ) {
                missingCpus += it.task.getResources().cpus
                tasks << it.task.taskId
            }
        }
        // cleanup the cache
        canRunTaskCache.clear()

        if( missingCpus && enabled ) {
            log.debug "### The following tasks have been waiting for more than $timeout -- required cpus=$missingCpus; taskIds=${tasks.join(',')}"
            requestNewCpus(missingCpus)
        }

    }

    /**
     * Request a specific amount of cpus spinning new instances in the cloud, making sure to not
     * overcome the current max instances limit
     *
     * @param cpus The number of cpus requested
     */
    @PackageScope
    void requestNewCpus( int cpus ) {
        def type = driver.describeInstanceType(scalerConfig.instanceType)
        if( !type ) {
            log.warn "### Can't find a instance type description: ${scalerConfig.instanceType}"
            return
        }

        def nodeNeeded = (int)Math.ceil( cpus / type.cpus )
        def currentSize = getClusterSize()
        def maxSize = scalerConfig.getMaxInstances()

        if( currentSize >= maxSize ) {
            log.debug "### Can't grow the cluster more, current size reached the cluster limit -- missing-cpus: $cpus; node-needed: $nodeNeeded; current-size: $currentSize; max-size: ${maxSize}"
            return
        }

        def num = Math.min( nodeNeeded, maxSize-currentSize )
        log.debug "### Requesting $num instance(s) of type: $type -- missing-cpus: $cpus; node-needed: $nodeNeeded; current-size: $currentSize; max-size: ${maxSize}"
        requestNewInstances(num)
    }

    /**
     * Submit the request for new cloud instances
     *
     * @param num The number of instances requested
     */
    @PackageScope
    void requestNewInstances( int num ) {
        // -- launch the requested nodes
        def ids = driver.launchInstances(num, scalerConfig)

        // -- disable cluster auto-scaling until all instance have joined the cluster
        this.enabled = false
        this.pendingInstanceIds = new ConcurrentLinkedQueue<>(ids)

        // -- tag the instance
        driver.waitInstanceStatus(ids, CloudInstanceStatus.STARTED)
        driver.tagInstances(ids, scalerConfig)

    }

    /**
     * Check if task resource requirements are fulfilled by a specified cloud instance type
     *
     * @param task A {@link IgBaseTask} instance modeling the task resources request
     * @param instanceType The instance type identifier e.g. {@code m4.xlarge}
     * @return {@code true} when the task can be execute in the specified instance type, {@code false} otherwise
     */
    @PackageScope
    boolean canRunOnNewInstance( IgBaseTask task, String instanceType ) {
        final taskId = task.getTaskId()
        final req = task.getResources()

        if( req.cpus == 1 && !req.memory )
            return true

        // -- make sure that tasks can be fulfilled by the new instance
        def total = driver.describeInstanceType(instanceType)
        if( !total ) {
            log.warn "### Unknown instance type: $instanceType"
            return false
        }

        if( req.cpus > total.cpus ) {
            log.warn("### Task (id=${taskId}) exceed the number of CPUs provided by autoscaling instance type: ${instanceType} -- req: ${req.cpus}; provided: ${total.cpus}")
            return false
        }
        if( req.memory && req.memory > total.memory ) {
            log.warn("### Task (id=${taskId}) exceed the amount of memory provided by autoscaling instance type: ${instanceType} -- req: ${req.memory}; provided: ${total.memory}")
            return false
        }

        return true
    }

    /**
     * Check a given task with some resources request can be executed in
     * any node in the current cluster topology
     *
     * @param task
     *      A {@link IgBaseTask} instance modeling the task requesting some computation resources
     * @return
     *      {@code 0} there's no free resources in the current cluster topology to fulfill the specified task
     *      {@code 1} there's at least one node with enough free resources to execute the task, and
     *      {@code 2} there's no node with enough to execute the task
     */
    @PackageScope
    byte canRunOnExistingNodes( IgBaseTask task ) {
        final taskId = task.getTaskId()
        final req = task.getResources()

        // keep the result in a local cache to avoid to make too many
        // queries over the distributed cache
        def key = (List)[ req.cpus, req.memory ]
        def found = canRunTaskCache.get(key)
        if( found != null )
            return found

        def overflow = true
        def fulfil = false
        try {
            def itr = workerNodes.values().iterator()
            while( itr.hasNext() ) {
                def node = itr.next()

                def total = node.resources
                if( req.cpus <= total.cpus && ( !req.memory || req.memory <= total.memory) ) {
                    // at least one instance has enough resources => reset the `overflow` flag
                    overflow = false
                }

                def free = node.free
                if( free && req.cpus <= free.cpus && (!req.memory || req.memory <= free.memory)) {
                    // at least one instance has enough free resources => set `fulfil` flag and exit
                    fulfil = true
                    break
                }
            }
        }
        catch ( Exception e ) {
            log.debug "### Oops.. Cannot establish resources availability: taskId=$taskId", e
        }

        if( overflow ) {
            log.warn "### Task (id=${taskId}) requests an amount of resources not available in any node in the current cluster topology -- CPUs: ${req.cpus}; memory: ${req.memory?:'-'}"
        }

        byte result = fulfil ? 1 : ( overflow ? 2 : 0 )
        canRunTaskCache.put(key, result)
        return result
    }

    /**
     * @return The ID of the local cluster node
     */
    @PackageScope getLocalNodeId() {
        ignite.cluster().localNode().id()
    }

    @PackageScope void checkIdleNodes() {

        if( !enabled || !scalerConfig.terminateWhenIdle )
            return

        final timeout = scalerConfig.getIdleTimeout()
        final idleNodeIds = []
        final itr = workerNodes.values().iterator()
        while( itr.hasNext() )  {
            final node = itr.next()

            try {
                if( node.isIdle(timeout) ) {
                    log.debug "### Idle node detected: $node"
                    idleNodeIds << node.nodeId
                }
            }
            catch( Exception e ) {
                log.debug "### Oops.. Failed to check idle node info: $node -- Cause: ${e.message ?: e}"
            }
        }

        // -- do not commit a suicide, remove this node
        def local = getLocalNodeId()
        idleNodeIds.remove(local)

        if( !idleNodeIds )
            return

        def killList = applyTerminationPolicy(idleNodeIds)
        if( killList ) {
            killNodes(killList)
        }
    }


    @PackageScope List<String> applyTerminationPolicy( Collection<UUID> idleNodeIds) {

        def killList = new ArrayList<String>(idleNodeIds.size())
        def itr = idleNodeIds.iterator()
        while( itr.hasNext() ) {

            final nodeId = itr.next()
            final data = workerNodes.get(nodeId)
            if( !data.instanceId ) {
                log.debug "### Oops.. Missing cloud instance id for node: [${data.hostName}] $nodeId -- Ignore termination request"
                continue
            }

            killList << data.instanceId
        }

       return killList
    }

    @PackageScope void killNodes(List<String> killList) {
        // avoid to terminate too many instances
        final instanceIds = ensureMinSize(killList)
        log.debug "### Killing instances: ids=${instanceIds.join(', ')}"
        requestTerminateInstances(instanceIds)
    }

    @PackageScope void requestTerminateInstances( List<String> instanceIds ) {
        driver.terminateInstances(instanceIds)
        driver.waitInstanceStatus(instanceIds, CloudInstanceStatus.TERMINATED)
    }

    @PackageScope
    List<String> ensureMinSize( List<String> killList ) {
        if( killList.size() >= clusterSize )
            throw new IllegalStateException("Can't kill all cluster nodes")

        int newSize = clusterSize - killList.size()
        def delta = newSize - scalerConfig.minInstances
        if( delta<0 ) {
            int last = delta-1
            return killList[0..last]
        }
        return killList
    }

    /**
     * Verifies periodically the status of pending instances requested
     * by the auto-scaler
     */
    @PackageScope void checkPendingInstances() {
        //TODO
    }

    /**
     * Shutdown the auto-scaler thread
     */
    @Override
    void close() throws IOException {
        watchdog.shutdownNow()
    }
}
