/*
 * Copyright (c) 2013-2015, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2015, Paolo Di Tommaso and the respective authors.
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
package nextflow.scheduler;

import static nextflow.scheduler.JobSchedulerSpi.THIEF_NODE_ATTR;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.UUID;

import org.apache.ignite.IgniteException;
import org.apache.ignite.IgniteLogger;
import org.apache.ignite.cluster.ClusterNode;
import org.apache.ignite.internal.util.typedef.F;
import org.apache.ignite.internal.util.typedef.internal.S;
import org.apache.ignite.internal.util.typedef.internal.U;
import org.apache.ignite.resources.LoggerResource;
import org.apache.ignite.spi.IgniteSpiAdapter;
import org.apache.ignite.spi.IgniteSpiConfiguration;
import org.apache.ignite.spi.IgniteSpiConsistencyChecked;
import org.apache.ignite.spi.IgniteSpiException;
import org.apache.ignite.spi.IgniteSpiMultipleInstancesSupport;
import org.apache.ignite.spi.failover.FailoverContext;
import org.apache.ignite.spi.failover.FailoverSpi;

/**
 * Implements the jobs failover strategy used by the nextflow scheduler.
 * This is basically a clone of the standard {@link org.apache.ignite.spi.failover.jobstealing.JobStealingFailoverSpi}
 * with the exception that patch the issue https://issues.apache.org/jira/browse/IGNITE-1267
 * allowing nodes not in the job topology to steal it.
 */
@IgniteSpiMultipleInstancesSupport(true)
@IgniteSpiConsistencyChecked(optional = true)
public class JobFailoverSpi extends IgniteSpiAdapter implements FailoverSpi,
        JobFailoverSpiMBean {
    /** Maximum number of attempts to execute a failed job on another node (default is {@code 5}). */
    public static final int DFLT_MAX_FAILOVER_ATTEMPTS = 5;

    /**
     * Name of job context attribute containing all nodes a job failed on. Note
     * that this list does not include nodes that a job was stolen from.
     *
     * @see org.apache.ignite.compute.ComputeJobContext
     */
    static final String FAILED_NODE_LIST_ATTR = "gg:failover:failednodelist";

    /**
     * Name of job context attribute containing current failover attempt count.
     * This count is incremented every time the same job gets failed over to
     * another node for execution if it was not successfully stolen.
     *
     * @see org.apache.ignite.compute.ComputeJobContext
     */
    static final String FAILOVER_ATTEMPT_COUNT_ATTR = "gg:failover:attemptcount";

    /** Maximum failover attempts job context attribute name. */
    private static final String MAX_FAILOVER_ATTEMPT_ATTR = "gg:failover:maxattempts";

    /** Injected grid logger. */
    @LoggerResource
    private IgniteLogger log;

    /** Maximum number of attempts to execute a failed job on another node. */
    private int maxFailoverAttempts = DFLT_MAX_FAILOVER_ATTEMPTS;

    /** Number of jobs that were failed over. */
    private int totalFailedOverJobs;

    /** Number of jobs that were stolen. */
    private int totalStolenJobs;

    /** {@inheritDoc} */
    @Override public int getMaximumFailoverAttempts() {
        return maxFailoverAttempts;
    }

    /**
     * Sets maximum number of attempts to execute a failed job on another node.
     * If job gets stolen and thief node exists then it is not considered as
     * failed job.
     * If not specified, {@link #DFLT_MAX_FAILOVER_ATTEMPTS} value will be used.
     * <p>
     * Note this value must be identical for all grid nodes in the grid.
     *
     * @param maxFailoverAttempts Maximum number of attempts to execute a failed
     *      job on another node.
     */
    @IgniteSpiConfiguration(optional = true)
    public void setMaximumFailoverAttempts(int maxFailoverAttempts) {
        this.maxFailoverAttempts = maxFailoverAttempts;
    }

    /** {@inheritDoc} */
    @Override public int getTotalFailedOverJobsCount() {
        return totalFailedOverJobs;
    }

    /** {@inheritDoc} */
    @Override public int getTotalStolenJobsCount() {
        return totalStolenJobs;
    }

    /** {@inheritDoc} */
    @Override public Map<String, Object> getNodeAttributes() throws IgniteSpiException {
        return F.<String, Object>asMap(createSpiAttributeName(MAX_FAILOVER_ATTEMPT_ATTR), maxFailoverAttempts);
    }

    /** {@inheritDoc} */
    @Override public void spiStart(String gridName) throws IgniteSpiException {
        // Start SPI start stopwatch.
        startStopwatch();

        assertParameter(maxFailoverAttempts >= 0, "maximumFailoverAttempts >= 0");

        if (log.isDebugEnabled())
            log.debug(configInfo("maxFailoverAttempts", maxFailoverAttempts));

        registerMBean(gridName, this, JobFailoverSpiMBean.class);

        // Ack ok start.
        if (log.isDebugEnabled())
            log.debug(startInfo());
    }

    /** {@inheritDoc} */
    @Override public void spiStop() throws IgniteSpiException {
        unregisterMBean();

        // Ack ok stop.
        if (log.isDebugEnabled())
            log.debug(stopInfo());
    }

    /** {@inheritDoc} */
    @SuppressWarnings("unchecked")
    @Override public ClusterNode failover(FailoverContext ctx, List<ClusterNode> top) {
        assert ctx != null;
        assert top != null;

        if (top.isEmpty()) {
            U.warn(log, "Received empty subgrid and is forced to fail.");

            // Nowhere to failover to.
            return null;
        }

        Integer failoverCnt = ctx.getJobResult().getJobContext().getAttribute(FAILOVER_ATTEMPT_COUNT_ATTR);

        if (failoverCnt == null)
            failoverCnt = 0;

        if (failoverCnt > maxFailoverAttempts) {
            U.error(log, "Failover count exceeded maximum failover attempts parameter [failedJob=" +
                    ctx.getJobResult().getJob() + ", maxFailoverAttempts=" + maxFailoverAttempts + ']');

            return null;
        }

        if (failoverCnt == maxFailoverAttempts) {
            U.warn(log, "Job failover failed because number of maximum failover attempts is exceeded [failedJob=" +
                    ctx.getJobResult().getJob() + ", maxFailoverAttempts=" + maxFailoverAttempts + ']');

            return null;
        }

        try {
            ClusterNode thief = null;
            boolean isNodeFailed = false;

            UUID thiefId = ctx.getJobResult().getJobContext().getAttribute(THIEF_NODE_ATTR);

            if (thiefId != null) {
                // Clear attribute.
                ctx.getJobResult().getJobContext().setAttribute(THIEF_NODE_ATTR, null);

                thief = getSpiContext().node(thiefId);

                if (thief != null) {
                    // If sender != receiver.
                    if (thief.equals(ctx.getJobResult().getNode())) {
                        U.error(log, "Job stealer node is equal to job node (will fail-over using load-balancing): " + thief.id());

                        isNodeFailed = true;

                        thief = null;
                    }
// -- disable this check due issue https://issues.apache.org/jira/browse/IGNITE-1267
//                    else if (!top.contains(thief)) {
//                        U.warn(log, "Thief node is not part of task topology  (will fail-over using load-balancing) " +
//                                "[thief=" + thiefId + ", topSize=" + top.size() + ']');
//
//                        thief = null;
//                    }

                    if (log.isTraceEnabled())
                        log.trace("Failing-over stolen job [from=" + ctx.getJobResult().getNode() + ", to=" + thief + ']');
                }
                else {
                    isNodeFailed = true;

                    U.warn(log, "Thief node left grid (will fail-over using load balancing): " + thiefId);
                }
            }
            else
                isNodeFailed = true;

            // If job was not stolen or stolen node is not part of topology,
            // then failover the regular way.
            if (thief == null) {
                Collection<UUID> failedNodes = ctx.getJobResult().getJobContext().getAttribute(FAILED_NODE_LIST_ATTR);

                if (failedNodes == null)
                    failedNodes = U.newHashSet(1);

                if (isNodeFailed)
                    failedNodes.add(ctx.getJobResult().getNode().id());

                // Set updated failed node set into job context.
                ctx.getJobResult().getJobContext().setAttribute(FAILED_NODE_LIST_ATTR, failedNodes);

                // Copy.
                List<ClusterNode> newTop = new ArrayList<>(top.size());

                for (ClusterNode n : top) {
                    // Add non-failed nodes to topology.
                    if (!failedNodes.contains(n.id()))
                        newTop.add(n);
                }

                if (newTop.isEmpty()) {
                    U.warn(log, "Received topology with only nodes that job had failed on (forced to fail) [failedNodes=" + failedNodes + ']');

                    // Nowhere to failover to.
                    return null;
                }

                thief = ctx.getBalancedNode(newTop);

                if (thief == null)
                    U.warn(log, "Load balancer returned null node for topology: " + newTop);
            }

            if (isNodeFailed)
                // This is a failover, not stealing.
                failoverCnt++;

            // Even if it was stealing and thief node left grid we assume
            // that it is failover because of the fail.
            ctx.getJobResult().getJobContext().setAttribute(FAILOVER_ATTEMPT_COUNT_ATTR, failoverCnt);

            if (thief != null) {
                totalFailedOverJobs++;

                if (isNodeFailed) {
                    U.warn(log, "Failed over job to a new node [newNode=" + thief.id() +
                            ", oldNode=" + ctx.getJobResult().getNode().id() +
                            ", sesId=" + ctx.getTaskSession().getId() +
                            ", job=" + ctx.getJobResult().getJob() +
                            ", jobCtx=" + ctx.getJobResult().getJobContext() +
                            ", task=" + ctx.getTaskSession().getTaskName() + ']');
                }
                else {
                    totalStolenJobs++;
                    if (log.isDebugEnabled())
                        log.debug("Stealing job to a new node [newNode=" + thief.id() +
                                ", oldNode=" + ctx.getJobResult().getNode().id() +
                                ", sesId=" + ctx.getTaskSession().getId() +
                                ", job=" + ctx.getJobResult().getJob() +
                                ", jobCtx=" + ctx.getJobResult().getJobContext() +
                                ", task=" + ctx.getTaskSession().getTaskName() + ']');
                }
            }

            return thief;
        }
        catch (IgniteException e) {
            U.error(log, "Failed to get next balanced node for failover: " + ctx, e);

            return null;
        }
    }

    /** {@inheritDoc} */
    @Override protected List<String> getConsistentAttributeNames() {
        return Collections.singletonList(createSpiAttributeName(MAX_FAILOVER_ATTEMPT_ATTR));
    }

    /** {@inheritDoc} */
    @Override public String toString() {
        return S.toString(JobFailoverSpi.class, this);
    }
}