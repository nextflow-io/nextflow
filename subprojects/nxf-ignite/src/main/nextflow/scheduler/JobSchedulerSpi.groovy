/*
 * Copyright (c) 2013-2016, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2016, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
 */


/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */


package nextflow.scheduler
import static org.apache.ignite.events.EventType.EVT_NODE_FAILED
import static org.apache.ignite.events.EventType.EVT_NODE_JOINED
import static org.apache.ignite.events.EventType.EVT_NODE_LEFT

import java.lang.reflect.Field
import java.util.concurrent.ConcurrentMap
import java.util.concurrent.atomic.AtomicInteger

import nextflow.daemon.IgGridFactory
import nextflow.executor.IgBaseTask
import nextflow.util.Duration
import org.apache.ignite.cluster.ClusterNode
import org.apache.ignite.compute.ComputeJobContext
import org.apache.ignite.events.DiscoveryEvent
import org.apache.ignite.events.Event
import org.apache.ignite.internal.managers.communication.GridMessageListener
import org.apache.ignite.internal.managers.eventstorage.GridLocalEventListener
import org.apache.ignite.internal.util.tostring.GridToStringInclude
import org.apache.ignite.internal.util.typedef.F
import org.apache.ignite.internal.util.typedef.internal.A
import org.apache.ignite.internal.util.typedef.internal.S
import org.apache.ignite.internal.util.typedef.internal.U
import org.apache.ignite.spi.IgniteSpiAdapter
import org.apache.ignite.spi.IgniteSpiConfiguration
import org.apache.ignite.spi.IgniteSpiConsistencyChecked
import org.apache.ignite.spi.IgniteSpiContext
import org.apache.ignite.spi.IgniteSpiException
import org.apache.ignite.spi.IgniteSpiMultipleInstancesSupport
import org.apache.ignite.spi.collision.CollisionContext
import org.apache.ignite.spi.collision.CollisionExternalListener
import org.apache.ignite.spi.collision.CollisionJobContext
import org.apache.ignite.spi.collision.CollisionSpi
import org.apache.ignite.spi.collision.jobstealing.JobStealingDisabled
import org.apache.ignite.spi.collision.jobstealing.JobStealingRequest
import org.jsr166.ConcurrentHashMap8
import org.jsr166.ConcurrentLinkedDeque8
import org.slf4j.Logger
import org.slf4j.LoggerFactory
/**
 * Implements the scheduling strategy for nextflow job. This code is largely based
 * on standard {@link org.apache.ignite.spi.collision.jobstealing.JobStealingCollisionSpi}.
 *
 * The main difference is that it handlers properly resources request declared at process
 * level such as cpus, memory and disk storage requirements
 *
 * Stealing logic:
 * - The {@link #checkIdle(ResourceContext)} checks if there are pending jobs to be executed
 *  in the current node. If there aren't wil send a `steal message` using the {@link #JOB_STEALING_COMM_TOPIC} topic.
 *
 * - The {@link #checkBusy(ResourceContext)} check for pending steal request, and cancel a job in waiting status
 *  and mark it as stolen adding the `thief` node information in the job context attributes
 *
 * - The {@link JobFailoverSpi} picks all rejected jobs. A jo marked as stolen is assigned to the `thief` node
 *   as defined in the job context attribute
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@SuppressWarnings("SynchronizationOnLocalVariableOrMethodParameter")
@IgniteSpiMultipleInstancesSupport(true)
@IgniteSpiConsistencyChecked(optional = true)
public class JobSchedulerSpi extends IgniteSpiAdapter implements CollisionSpi, JobSchedulerSpiMBean {

    private static final Logger log = LoggerFactory.getLogger(JobSchedulerSpi.class);

    /** Maximum number of attempts to steal job by another node (default is {@code 5}). */
    public static final int DFLT_MAX_STEALING_ATTEMPTS = 5;

    public static final long DFLT_MAX_STEALING_EXPIRE_TIME = 5 * 60 * 1000 // 5 minutes

    /**
     * Default number of parallel jobs allowed (value is {@code 95} which is
     * slightly less same as default value of threads in the execution thread pool
     * to allow some extra threads for system processing).
     */
    public static final int DFLT_ACTIVE_JOBS_THRESHOLD = 95;

    /**
     * Default steal message expire time in milliseconds (value is {@code 1000}).
     * Once this time is elapsed and no response for steal message is received,
     * the message is considered lost and another steal message will be generated,
     * potentially to another node.
     */
    public static final long DFLT_MSG_EXPIRE_TIME = 1000;

    /**
     * Default threshold of waiting jobs. If number of waiting jobs exceeds this threshold,
     * then waiting jobs will become available to be stolen (value is {@code 0}).
     */
    public static final int DFLT_WAIT_JOBS_THRESHOLD = 0;

    /** Default start value for job priority (value is {@code 0}). */
    public static final int DFLT_JOB_PRIORITY = 0;

    /** Communication topic. */
    private static final String JOB_STEALING_COMM_TOPIC = "ignite.collision.job.stealing.topic";

    /** Job context attribute for storing thief node UUID (this attribute is used in job stealing failover SPI). */
    public static final String THIEF_NODE_ATTR = "ignite.collision.thief.node";

    /** Threshold of maximum jobs on waiting queue. */
    public static final String WAIT_JOBS_THRESHOLD_NODE_ATTR = "ignite.collision.wait.jobs.threshold";

    /** Threshold of maximum jobs executing concurrently. */
    public static final String ACTIVE_JOBS_THRESHOLD_NODE_ATTR = "ignite.collision.active.jobs.threshold";

    /**
     * Name of job context attribute containing current stealing attempt count.
     * This count is incremented every time the same job gets stolen for
     * execution.
     *
     * @see org.apache.ignite.compute.ComputeJobContext
     */
    public static final String STEALING_ATTEMPT_COUNT_ATTR = "ignite.stealing.attempt.count";

    /** Maximum stealing attempts attribute name. */
    public static final String MAX_STEALING_ATTEMPT_ATTR = "ignite.stealing.max.attempts";

    /** Stealing request expiration time attribute name. */
    public static final String MSG_EXPIRE_TIME_ATTR = "ignite.stealing.msg.expire.time";

    /** Stealing priority attribute name. */
    public static final String STEALING_PRIORITY_ATTR = "ignite.stealing.priority";

    /** The timestamp when a stealing was last attempt */
    public static final String STEALING_ATTEMPT_TIMESTAMP_ATTR = 'ignite.stealing.timestamp'

    /** Number of jobs that can be executed in parallel. */
    private volatile int activeJobsThreshold = DFLT_ACTIVE_JOBS_THRESHOLD;

    /** Configuration parameter defining waiting job count threshold for stealing to start. */
    @SuppressWarnings("RedundantFieldInitialization")
    private volatile int waitJobsThreshold = DFLT_WAIT_JOBS_THRESHOLD;

    /** Message expire time configuration parameter. */
    private volatile long msgExpireTime = DFLT_MSG_EXPIRE_TIME;

    /** Maximum number of attempts to steal job by another node. */
    private volatile int maxStealingAttempts = DFLT_MAX_STEALING_ATTEMPTS;

    private volatile long maxStealingExpireTime = DFLT_MAX_STEALING_EXPIRE_TIME

    /** Flag indicating whether job stealing is enabled. */
    private volatile boolean isStealingEnabled = true;

    /** Steal attributes. */
    @GridToStringInclude
    private Map<String, ? extends Serializable> stealAttrs;

    /** Number of jobs that were active last time. */
    private volatile int runningNum;

    /** Number of jobs that were waiting for execution last time. */
    private volatile int waitingNum;

    /** Number of currently held jobs. */
    private volatile int heldNum;

    /** Total number of stolen jobs. */
    private final AtomicInteger totalStolenJobsNum = new AtomicInteger();

    /** Map of sent messages. */
    private final ConcurrentMap<UUID, MessageInfo> sndMsgMap = new ConcurrentHashMap8<>();

    /** Map of received messages. */
    private final ConcurrentMap<UUID, MessageInfo> rcvMsgMap = new ConcurrentHashMap8<>();

    /** */
    private final Queue<ClusterNode> nodeQueue = new ConcurrentLinkedDeque8<>();

    /** */
    private CollisionExternalListener extLsnr;

    /** Discovery listener. */
    private GridLocalEventListener discoLsnr;

    /** Communication listener. */
    private GridMessageListener msgLsnr;

    /** Number of steal requests. */
    private final AtomicInteger stealReqs = new AtomicInteger();

    /** */
    private Comparator<CollisionJobContext> cmp;

    /** {@inheritDoc} */
    @IgniteSpiConfiguration(optional = true)
    @Override public void setActiveJobsThreshold(int activeJobsThreshold) {
        A.ensure(activeJobsThreshold >= 0, "activeJobsThreshold >= 0");

        this.activeJobsThreshold = activeJobsThreshold;
    }

    /** {@inheritDoc} */
    @Override public int getActiveJobsThreshold() {
        return activeJobsThreshold;
    }

    /** {@inheritDoc} */
    @IgniteSpiConfiguration(optional = true)
    @Override public void setWaitJobsThreshold(int waitJobsThreshold) {
        A.ensure(waitJobsThreshold >= 0, "waitJobsThreshold >= 0");

        this.waitJobsThreshold = waitJobsThreshold;
    }

    /** {@inheritDoc} */
    @Override public int getWaitJobsThreshold() {
        return waitJobsThreshold;
    }

    /** {@inheritDoc} */
    @IgniteSpiConfiguration(optional = true)
    @Override public void setMessageExpireTime(long msgExpireTime) {
        A.ensure(msgExpireTime > 0, "messageExpireTime > 0");

        this.msgExpireTime = msgExpireTime;
    }

    /** {@inheritDoc} */
    @Override public long getMessageExpireTime() {
        return msgExpireTime;
    }

    /** {@inheritDoc} */
    @IgniteSpiConfiguration(optional = true)
    @Override public void setStealingEnabled(boolean isStealingEnabled) {
        this.isStealingEnabled = isStealingEnabled;
    }

    /** {@inheritDoc} */
    @Override public boolean isStealingEnabled() {
        return isStealingEnabled;
    }

    /** {@inheritDoc} */
    @IgniteSpiConfiguration(optional = true)
    @Override public void setMaximumStealingAttempts(int maxStealingAttempts) {
        A.ensure(maxStealingAttempts > 0, "maxStealingAttempts > 0");

        this.maxStealingAttempts = maxStealingAttempts;
    }

    /** {@inheritDoc} */
    @Override public int getMaximumStealingAttempts() {
        return maxStealingAttempts;
    }

    public void setMaximumStealingExpireTime(Duration time ) {
        this.maxStealingExpireTime = time ? time.toMillis() : DFLT_MAX_STEALING_EXPIRE_TIME
    }

    public void getMaximumStealingExpireTime() {
        this.maxStealingExpireTime
    }


    /**
     * Configuration parameter to enable stealing to/from only nodes that
     * have these attributes set (see {@link org.apache.ignite.cluster.ClusterNode#attribute(String)} and
     * {@link org.apache.ignite.configuration.IgniteConfiguration#getUserAttributes()} methods).
     *
     * @param stealAttrs Node attributes to enable job stealing for.
     */
    @IgniteSpiConfiguration(optional = true)
    public void setStealingAttributes(Map<String, ? extends Serializable> stealAttrs) {
        this.stealAttrs = stealAttrs;
    }

    /** {@inheritDoc} */
    @Override public Map<String, ? extends Serializable> getStealingAttributes() {
        return stealAttrs;
    }

    /** {@inheritDoc} */
    @Override public int getCurrentRunningJobsNumber() {
        return runningNum;
    }

    /** {@inheritDoc} */
    @Override public int getCurrentHeldJobsNumber() {
        return heldNum;
    }

    /** {@inheritDoc} */
    @Override public int getCurrentWaitJobsNumber() {
        return waitingNum;
    }

    /** {@inheritDoc} */
    @Override public int getCurrentActiveJobsNumber() {
        return runningNum + heldNum;
    }

    /** {@inheritDoc} */
    @Override public int getTotalStolenJobsNumber() {
        return totalStolenJobsNum.get();
    }

    /** {@inheritDoc} */
    @Override public int getCurrentJobsToStealNumber() {
        return stealReqs.get();
    }

    /** {@inheritDoc} */
    @Override public Map<String, Object> getNodeAttributes() throws IgniteSpiException {
        return F.<String, Object>asMap(
                createSpiAttributeName(WAIT_JOBS_THRESHOLD_NODE_ATTR), waitJobsThreshold,
                createSpiAttributeName(ACTIVE_JOBS_THRESHOLD_NODE_ATTR), activeJobsThreshold,
                createSpiAttributeName(MAX_STEALING_ATTEMPT_ATTR), maxStealingAttempts,
                createSpiAttributeName(MSG_EXPIRE_TIME_ATTR), msgExpireTime);
    }

    /** {@inheritDoc} */
    @Override public void spiStart(String gridName) throws IgniteSpiException {
        assertParameter(activeJobsThreshold >= 0, "activeJobsThreshold >= 0");
        assertParameter(waitJobsThreshold >= 0, "waitJobsThreshold >= 0");
        assertParameter(msgExpireTime > 0, "messageExpireTime > 0");
        assertParameter(maxStealingAttempts > 0, "maxStealingAttempts > 0");

        // Start SPI start stopwatch.
        startStopwatch();

        // Ack parameters.
        //log.debug(configInfo("activeJobsThreshold", activeJobsThreshold));
        log.debug1(configInfo("waitJobsThreshold", waitJobsThreshold));
        log.debug1(configInfo("messageExpireTime", msgExpireTime));
        log.debug1(configInfo("maxStealingAttempts", maxStealingAttempts));
        log.debug1(configInfo("maxStealingExpireTime", maxStealingExpireTime));

        registerMBean(gridName, this, JobSchedulerSpiMBean.class);

        // Ack start.
        log.trace1(startInfo());
    }

    /** {@inheritDoc} */
    @Override public void spiStop() throws IgniteSpiException {
        unregisterMBean();

        // Ack ok stop.
        log.trace1(stopInfo());
    }

    /** {@inheritDoc} */
    @Override public void setExternalCollisionListener(CollisionExternalListener extLsnr) {
        this.extLsnr = extLsnr;
    }

    /** {@inheritDoc} */
    @Override protected void onContextInitialized0(IgniteSpiContext spiCtx) throws IgniteSpiException {
        spiCtx.addLocalEventListener(
                discoLsnr = new GridLocalEventListener() {
                    @SuppressWarnings("fallthrough")
                    @Override public void onEvent(Event evt) {
                        assert evt instanceof DiscoveryEvent;

                        DiscoveryEvent discoEvt = (DiscoveryEvent)evt;

                        UUID evtNodeId = discoEvt.eventNode().id();

                        switch (discoEvt.type()) {
                            case EVT_NODE_JOINED:
                                ClusterNode node = getSpiContext().node(evtNodeId);

                                if (node != null) {
                                    log.debug1("Node joined > id: " + evtNodeId );
                                    nodeQueue.offer(node);

                                    sndMsgMap.putIfAbsent(node.id(), new MessageInfo());
                                    rcvMsgMap.putIfAbsent(node.id(), new MessageInfo());
                                }

                                break;

                            case EVT_NODE_LEFT:
                            case EVT_NODE_FAILED:
                                log.debug1("Node left > id: " + evtNodeId );
                                Iterator<ClusterNode> iter = nodeQueue.iterator();

                                while (iter.hasNext()) {
                                    ClusterNode nextNode = iter.next();

                                    if (nextNode.id().equals(evtNodeId))
                                        iter.remove();
                                }

                                sndMsgMap.remove(evtNodeId);
                                rcvMsgMap.remove(evtNodeId);

                                break;

                            default:
                                assert false : "Unexpected event: " + evt;
                        }
                    }
                },
                EVT_NODE_FAILED,
                EVT_NODE_JOINED,
                EVT_NODE_LEFT
        );

        Collection<ClusterNode> rmtNodes = spiCtx.remoteNodes();

        for (ClusterNode node : rmtNodes) {
            UUID id = node.id();

            if (spiCtx.node(id) != null) {
                sndMsgMap.putIfAbsent(id, new MessageInfo());
                rcvMsgMap.putIfAbsent(id, new MessageInfo());

                // Check if node has concurrently left.
                if (spiCtx.node(id) == null) {
                    sndMsgMap.remove(id);
                    rcvMsgMap.remove(id);
                }
            }
        }

        nodeQueue.addAll(rmtNodes);

        Iterator<ClusterNode> iter = nodeQueue.iterator();

        while (iter.hasNext()) {
            ClusterNode nextNode = iter.next();

            if (spiCtx.node(nextNode.id()) == null)
                iter.remove();
        }

        /*
         * Receive a job steal request and add it to the `stealReqs` counter
         */
        spiCtx.addMessageListener(
                msgLsnr = new GridMessageListener() {
                    @Override public void onMessage(UUID nodeId, Object msg) {
                        log.trace1("Receiving stealing message: " + msg);

                        MessageInfo info = rcvMsgMap.get(nodeId);

                        if (info == null) {
                            log.trace1("Ignoring message steal request as discovery event has not yet been received for node: " + nodeId);
                            return;
                        }

                        int stealReqs0;

                        synchronized (info) {
                            JobStealingRequest req = (JobStealingRequest)msg;

                            // Increment total number of steal requests.
                            // Note that it is critical to increment total
                            // number of steal requests before resetting message info.
                            int delta = requestDelta(req);
                            stealReqs0 = stealReqs.addAndGet(delta - info.jobsToSteal());

                            info.reset(delta);
                        }

                        log.trace1("Received steal request [nodeId=${nodeId}, msg=${msg}, stealReqs=${stealReqs0}]");

                        CollisionExternalListener tmp = extLsnr;

                        // Let grid know that collisions should be resolved.
                        if (tmp != null)
                            tmp.onExternalCollision();
                    }
                },
                JOB_STEALING_COMM_TOPIC);
    }

    /** {@inheritDoc} */
    @Override public void onContextDestroyed0() {
        if (discoLsnr != null)
            getSpiContext().removeLocalEventListener(discoLsnr);

        if (msgLsnr != null)
            getSpiContext().removeMessageListener(msgLsnr, JOB_STEALING_COMM_TOPIC);
    }

    /** {@inheritDoc} */
    @Override
    public synchronized void onCollision(CollisionContext ctx) {
        assert ctx != null;

        ResourceContext res = new ResourceContext(ctx);
        heldNum = ctx.heldJobs().size();

        // save the current context
        saveNodeResources(res);
        activeJobsThreshold = res.getAvailCpus();

        // Check if there are any jobs to activate or reject.
        int rejected = checkBusy(res);

        totalStolenJobsNum.addAndGet(rejected);

        // No point of stealing jobs if some jobs were rejected.
        if (rejected > 0) {
            log.trace1("Total count of rejected jobs: " + rejected);
            return;
        }

        if (isStealingEnabled)
            // Check if there are jobs to steal.
            checkIdle(res);
    }

    protected void saveNodeResources(ResourceContext res) {
        final IgniteSpiContext ctx = getSpiContext();
        ctx.put(IgGridFactory.RESOURCE_CACHE, ctx.localNode().id(), res, 0);
    }

    protected ResourceContext getNodeResources(UUID nodeId) {
        return getSpiContext().get(IgGridFactory.RESOURCE_CACHE, nodeId);
    }

    /**
     * Check if node is busy and activate/reject proper number of jobs.
     *
     * @param res Resource.
     * @return Number of rejected jobs.
     */
    private int checkBusy(ResourceContext res) {

        int activeSize = res.getActiveJobs().size();
        int waitSize = res.getWaitingJobs().size();
        log.debug1("Scheduler jobs > active: ${activeSize} - pending: ${waitSize} - steal-reqs: ${stealReqs.get()}");

        waitingNum = res.getWaitingJobs().size();
        runningNum = activeSize;

        IgniteSpiContext ctx = getSpiContext();

        int activated = 0;
        int rejected = 0;

        Collection<CollisionJobContext> waitPriJobs = sortJobs(res.getWaitingJobs(), waitSize);

        for (CollisionJobContext waitCtx : waitPriJobs) {

            /*
             * try to activate a job
             */
            if (res.canActivate(waitCtx)) {
                // If job was activated/cancelled by another thread, then
                // this method is no-op.
                // We also need to make sure that job is not being rejected by another thread.
                if( waitCtx.activate() ) {
                    res.consumeTaskResource(waitCtx);
                    activated = res.getBusyCpus();
                    saveNodeResources(res);
                }
            }

            /*
             * check for stealing requests and offer jobs by rejecting (cancel) them
             */
            else {
                log.debug1( "Can't activate pending task -- Cause: ${res.message}" )

                if (stealReqs.get() > 0) {
                    if (waitCtx.getJob().getClass().isAnnotationPresent(JobStealingDisabled.class))
                        continue;

                    // Collision count attribute.
                    Integer stealingCnt = waitCtx.getJobContext().getAttribute(STEALING_ATTEMPT_COUNT_ATTR);
                    Long stealingTimestamp = waitCtx.getJobContext().getAttribute(STEALING_ATTEMPT_TIMESTAMP_ATTR)

                    // Check that maximum stealing attempt threshold
                    // has not been exceeded.
                    if (stealingCnt != null) {
                        // If job exceeded failover threshold, skip it.
                        if (stealingCnt >= maxStealingAttempts) {
                            long delta = stealingTimestamp ? System.currentTimeMillis()-stealingTimestamp : 0
                            if( delta < maxStealingExpireTime ) {
                                log.trace1("Waiting job exceeded stealing attempts and won't be rejected (will try other jobs on waiting list): " + waitCtx);
                                continue;
                            }
                            else {
                                // reset the stealing count
                                log.trace1("Resetting stealing count for job: ${waitCtx}")
                                waitCtx.getJobContext().setAttribute(STEALING_ATTEMPT_COUNT_ATTR, 0);
                                waitCtx.getJobContext().setAttribute(STEALING_ATTEMPT_TIMESTAMP_ATTR, null)
                                stealingCnt = 0;
                            }
                        }
                    }
                    else
                        stealingCnt = 0;

                    // Check if allowed to reject job.
                    int jobsToReject = waitPriJobs.size() - activated - rejected - waitJobsThreshold;

                    log.trace1("Jobs to reject count [jobsToReject=${jobsToReject}, waitCtx=${waitCtx}]");

                    if (jobsToReject <= 0)
                        break;

                    final Integer pri = waitCtx.getJobContext().getAttribute(STEALING_PRIORITY_ATTR) ?: DFLT_JOB_PRIORITY

                    // If we have an excess of waiting jobs, reject as many as there are
                    // requested to be stolen. Note, that we use lose total steal request
                    // counter to prevent excessive iteration over nodes under load.
                    for (Iterator<Map.Entry<UUID, MessageInfo>> iter = rcvMsgMap.entrySet().iterator();
                         iter.hasNext() && stealReqs.get() > 0;) {
                        Map.Entry<UUID, MessageInfo> entry = iter.next();

                        UUID nodeId = entry.getKey();

                        // Node has left topology.
                        if (ctx.node(nodeId) == null) {
                            iter.remove();
                            continue;
                        }

                        // give only jobs for which there's are enough
                        // computing resources in the thief node
                        ResourceContext thiefRes = getNodeResources(nodeId);
                        if( !thiefRes ) {
                            log.trace1("Decline steal request from node `$nodeId` due to missing resources context")
                            continue
                        }
                        if( !thiefRes.canActivate(waitCtx) ) {
                            log.trace1("Decline steal request from node `$nodeId` -- Cause: $res.message");
                            continue;
                        }

                        MessageInfo info = entry.getValue();

                        synchronized (info) {
                            int jobsAsked = info.jobsToSteal();
                            assert jobsAsked >= 0;

                            // Skip nodes that have not asked for jobs to steal.
                            if (jobsAsked == 0)
                            // Move to next node.
                                continue;

                            // If message is expired, ignore it.
                            if (info.expired()) {
                                // Subtract expired messages.
                                stealReqs.addAndGet(-info.jobsToSteal());
                                info.reset(0);
                                continue;
                            }

// -- disable this check due issue https://issues.apache.org/jira/browse/IGNITE-1267
//                        // Check that waiting job has thief node in topology.
//                        boolean found = false;
//
//                        for (UUID id : waitCtx.getTaskSession().getTopology()) {
//                            if (id.equals(nodeId)) {
//                                found = true;
//
//                                break;
//                            }
//                        }
//
//                        if (!found) {
//                            if (log.isDebugEnabled())
//                                log.debug("Thief node does not belong to task topology [thief=" + nodeId +
//                                        ", task=" + waitCtx.getTaskSession() + ']');
//
//                            continue;
//                        }

                            if (stealReqs.get() <= 0)
                                break;

                            // Need to make sure that job is not being
                            // rejected by another thread.
                            boolean cancel = waitCtx.getJobContext().getAttribute(THIEF_NODE_ATTR) == null;

                            if (cancel) {
                                // Mark job as stolen.
                                waitCtx.getJobContext().setAttribute(THIEF_NODE_ATTR, nodeId);
                                waitCtx.getJobContext().setAttribute(STEALING_ATTEMPT_COUNT_ATTR, stealingCnt + 1);
                                waitCtx.getJobContext().setAttribute(STEALING_PRIORITY_ATTR, pri + 1);
                                waitCtx.getJobContext().setAttribute(STEALING_ATTEMPT_TIMESTAMP_ATTR, System.currentTimeMillis())

                                log.trace1("Will try to reject job due to steal request [ctx=${waitCtx}, thief=${nodeId}]");

                                int i = stealReqs.decrementAndGet();

                                // reject the job because there's a steal request
                                // the `JobFailOverSpi` will pick it assign to the `thief` node
                                if (i >= 0 && waitCtx.cancel()) {
                                    rejected++;

                                    info.reset(jobsAsked - 1);

                                    log.trace1("Rejected job due to steal request [ctx=${waitCtx}, nodeId=${nodeId}]");
                                }
                                else {
                                    log.trace1("Failed to reject job [i=${i}]");

                                    waitCtx.getJobContext().setAttribute(THIEF_NODE_ATTR, null);
                                    waitCtx.getJobContext().setAttribute(STEALING_ATTEMPT_COUNT_ATTR, stealingCnt);
                                    waitCtx.getJobContext().setAttribute(STEALING_PRIORITY_ATTR, pri);
                                    waitCtx.getJobContext().setAttribute(STEALING_ATTEMPT_TIMESTAMP_ATTR, null)

                                    stealReqs.incrementAndGet();
                                }
                            }

                            // Move to next job.
                            break;
                        }
                    }
                }
                else
                // No more jobs to steal or activate.
                    break;
            }

        }

        return rejected;
    }

    /**
     * Sort jobs by priority from high to lowest value.
     *
     * @param waitJobs Waiting jobs.
     * @param waitSize Snapshot size.
     * @return Sorted waiting jobs by priority.
     */
    private Collection<CollisionJobContext> sortJobs(Collection<CollisionJobContext> waitJobs, int waitSize) {
        List<CollisionJobContext> passiveList = new ArrayList<>(waitJobs.size());

        int i = 0;

        for (CollisionJobContext waitJob : waitJobs) {
            passiveList.add(waitJob);

            if (i++ == waitSize)
                break;
        }

        Collections.sort(passiveList, comparator());

        return passiveList;
    }

    /**
     * @return Comparator.
     */
    private Comparator<CollisionJobContext> comparator() {
        if (cmp == null) {
            cmp = new Comparator<CollisionJobContext>() {
                @Override public int compare(CollisionJobContext o1, CollisionJobContext o2) {
                    int p1 = getJobPriority(o1.getJobContext());
                    int p2 = getJobPriority(o2.getJobContext());

                    return Integer.compare(p2, p1);
                }
            };
        }

        return cmp;
    }

    /**
     * Gets job priority from task context. If job has no priority default one will be used.
     *
     * @param ctx Job context.
     * @return Job priority.
     */
    static private int getJobPriority(ComputeJobContext ctx) {
        assert ctx != null;

        Integer p;

        try {
            p = ctx.getAttribute(STEALING_PRIORITY_ATTR);
        }
        catch (ClassCastException e) {
            log.error("Type of job context priority attribute '" + STEALING_PRIORITY_ATTR +
                    "' is not java.lang.Integer (will use default priority) [type=" +
                    ctx.getAttribute(STEALING_PRIORITY_ATTR).getClass() + ", dfltPriority=" + DFLT_JOB_PRIORITY + ']', e);

            p = DFLT_JOB_PRIORITY;
        }

        if (p == null)
            p = DFLT_JOB_PRIORITY;

        return p;
    }

    /**
     * Check if the node is idle and steal as many jobs from other nodes
     * as possible.
     *
     * @param res Active jobs.
     */
    private void checkIdle(ResourceContext res) {
        // Check for overflow.
        int max = waitJobsThreshold + activeJobsThreshold;

        if (max < 0)
            max = Integer.MAX_VALUE;

        Collection<CollisionJobContext> waitJobs = res.getWaitingJobs();
        Collection<CollisionJobContext> activeJobs = res.getActiveJobs();
        int jobsToSteal = max - (waitJobs.size() + activeJobs.size());

        log.trace1("Total number of jobs to be stolen: " + jobsToSteal);

        if (jobsToSteal > 0) {
            int jobsLeft = jobsToSteal;

            ClusterNode next;

            int nodeCnt = getSpiContext().remoteNodes().size();

            int idx = 0;

            while (jobsLeft > 0 && idx++ < nodeCnt && (next = nodeQueue.poll()) != null) {
                if (getSpiContext().node(next.id()) == null)
                    continue;

                // Remote node does not have attributes - do not steal from it.
                if (stealAttrs && (next.attributes() == null || !U.containsAll(next.attributes(), stealAttrs))) {
                    log.trace1("Skip node as it does not have all attributes: " + next.id());

                    continue;
                }

                int delta = 0;

                try {
                    MessageInfo msgInfo = sndMsgMap.get(next.id());

                    if (msgInfo == null) {
                        log.trace1("Failed to find message info for node: " + next.id());
                        // Node left topology or SPI has not received message for it.
                        continue;
                    }

                    Integer waitThreshold =
                            next.attribute(createSpiAttributeName(WAIT_JOBS_THRESHOLD_NODE_ATTR));

                    if (waitThreshold == null) {
                        log.error("Remote node is not configured with GridJobStealingCollisionSpi and " +
                                "jobs will not be stolen from it (you must stop it and update its configuration to use " +
                                "GridJobStealingCollisionSpi): " + next)
                        continue;
                    }

                    delta = next.metrics().getCurrentWaitingJobs() - waitThreshold;

                    log.trace1("Maximum number of jobs to steal from node [jobsToSteal=${delta}, node=${next.id()}]");

                    // Nothing to steal from this node.
                    if (delta <= 0)
                        continue;

                    synchronized (msgInfo) {
                        if (!msgInfo.expired() && msgInfo.jobsToSteal() > 0) {
                            // Count messages being waited for as present.
                            jobsLeft -= msgInfo.jobsToSteal();

                            continue;
                        }

                        if (jobsLeft < delta)
                            delta = jobsLeft;

                        jobsLeft -= delta;

                        msgInfo.reset(delta);
                    }

                    // Send request to remote node to steal jobs.
                    // Message is a plain integer represented by 'delta'.
                    log.trace1("Sending message to steal ${delta} jobs");
                    getSpiContext().send(next, newStealRequest(delta), JOB_STEALING_COMM_TOPIC);
                }
                catch (IgniteSpiException e) {
                    log.error("Failed to send job stealing message to node: " + next, e);

                    // Rollback.
                    jobsLeft += delta;
                }
                finally {
                    // If node is alive, add back to the end of the queue.
                    if (getSpiContext().node(next.id()) != null)
                        nodeQueue.offer(next);
                }
            }
        }
    }

    /** {@inheritDoc} */
    @Override protected List<String> getConsistentAttributeNames() {
        List<String> attrs = new ArrayList<>(2);

        attrs.add(createSpiAttributeName(MAX_STEALING_ATTEMPT_ATTR));
        attrs.add(createSpiAttributeName(MSG_EXPIRE_TIME_ATTR));

        return attrs;
    }

    /** {@inheritDoc} */
    @Override public String toString() {
        return S.toString(JobSchedulerSpi.class, this);
    }

    /**
     * Creates a new instance of {@link JobStealingRequest} setting the delta value
     * properly. This hack is required because the `delta` field is declared private.
     *
     * @param value The request delta value
     * @return A new instance of {@link JobStealingRequest}
     */
    private JobStealingRequest newStealRequest(int value) {
        JobStealingRequest req = new JobStealingRequest();
        try {
            Field delta = JobStealingRequest.class.getDeclaredField("delta");
            delta.setAccessible(true);
            delta.set(req, value);
            return req;
        }
        catch( Exception e ) {
            throw new IllegalStateException("Cannot create steal request", e);
        }
    }

    /**
     * This hack is required because the `delta` field is declared private.
     */
    private int requestDelta(JobStealingRequest req) {
        try {
            Field delta = JobStealingRequest.class.getDeclaredField("delta");
            delta.setAccessible(true);
            return (Integer)delta.get(req);
        }
        catch( Exception e ) {
            throw new IllegalStateException("Cannot access steal request delta", e);
        }
    }

    /**
     *
     */
    private class MessageInfo {
        /** */
        private int jobsToSteal;

        /** */
        private long ts = U.currentTimeMillis();

        /**
         * @return Job to steal.
         */
        int jobsToSteal() {
            assert Thread.holdsLock(this);

            return jobsToSteal;
        }

        /**
         * @return {@code True} if message is expired.
         */
        boolean expired() {
            assert Thread.holdsLock(this);

            return jobsToSteal > 0 && U.currentTimeMillis() - ts >= msgExpireTime;
        }

        /**
         * @param jobsToSteal Jobs to steal.
         */
        void reset(int jobsToSteal) {
            assert Thread.holdsLock(this);

            this.jobsToSteal = jobsToSteal;

            ts = U.currentTimeMillis();
        }

        /** {@inheritDoc} */
        @Override public String toString() {
            return S.toString(MessageInfo.class, this);
        }
    }

}