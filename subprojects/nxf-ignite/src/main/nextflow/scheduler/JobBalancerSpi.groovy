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
import groovy.transform.CompileStatic
import nextflow.daemon.IgGridFactory
import nextflow.executor.IgBaseTask
import org.apache.ignite.cluster.ClusterNode
import org.apache.ignite.compute.ComputeJob
import org.apache.ignite.compute.ComputeTaskSession
import org.apache.ignite.internal.util.typedef.internal.A
import org.apache.ignite.internal.util.typedef.internal.S
import org.apache.ignite.spi.IgniteSpiAdapter
import org.apache.ignite.spi.IgniteSpiConsistencyChecked
import org.apache.ignite.spi.IgniteSpiException
import org.apache.ignite.spi.IgniteSpiMultipleInstancesSupport
import org.apache.ignite.spi.loadbalancing.LoadBalancingSpi
import org.jetbrains.annotations.Nullable
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@IgniteSpiMultipleInstancesSupport(true)
@IgniteSpiConsistencyChecked(optional = true)
public class JobBalancerSpi extends IgniteSpiAdapter implements LoadBalancingSpi, JobBalancerSpiMBean {

    private static final Logger log = LoggerFactory.getLogger(JobBalancerSpi)

    /** Random number generator. */
    private static final Random RAND = new Random();

    /** {@inheritDoc} */
    @Override public void spiStart(@Nullable String gridName) throws IgniteSpiException {
        startStopwatch();

        registerMBean(gridName, this, JobBalancerSpiMBean.class);

        // Ack ok start.
        log.debug1(startInfo());
    }

    /** {@inheritDoc} */
    @Override public void spiStop() throws IgniteSpiException {
        unregisterMBean();

        // Ack ok stop.
        log.debug1(stopInfo());
    }


    /** {@inheritDoc} */
    @Override public ClusterNode getBalancedNode(ComputeTaskSession ses, List<ClusterNode> top, ComputeJob job) {
        A.notNull(ses, "ses");
        A.notNull(top, "top");
        A.notNull(job, "job");

        List<ClusterNode> copy = new ArrayList<>(top)
        while( copy.size() ) {
            final idx = RAND.nextInt(copy.size())
            ClusterNode node = copy.get(idx);
            copy.remove(idx)
            def res = (ResourceContext)getSpiContext().get(IgGridFactory.RESOURCE_CACHE, node.id())

            if( job instanceof IgBaseTask && res ) {
                log.trace1 "Balancer checking node: ${node.id()}  -- Resources provided: $res -- Resources requested: ${job.resources}"

                if( job.resources.cpus > res.freeCpus ) continue
                if( job.resources.memory && job.resources.memory > res.freeMemory ) continue
                if( job.resources.disk && job.resources.disk > res.freeDisk ) continue

                log.trace1 "Balancer picked node: ${node.id()}"
                return node
            }

            return node
        }

        ClusterNode node = top.get(RAND.nextInt(top.size()))
        log.trace1 "Balancer cannot determine a balanced node -- A random one was picked > $node"
        return node
    }


    /** {@inheritDoc} */
    @Override public String toString() {
        return S.toString(JobBalancerSpi.class, this);
    }
}