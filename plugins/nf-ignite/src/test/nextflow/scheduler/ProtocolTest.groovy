/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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
import nextflow.executor.IgBaseTask
import nextflow.scheduler.Protocol.Resources
import nextflow.scheduler.Protocol.TaskHolder
import nextflow.util.ClusterConfig
import nextflow.util.Duration
import nextflow.util.MemoryUnit
import nextflow.util.SysHelper
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ProtocolTest extends Specification {

    def 'should clone resources obj' () {

        def res = new Resources(cpus:4, memory: MemoryUnit.of('8G'), disk: MemoryUnit.of('100GB'))
        when:
        def copy = res.clone()
        then:
        res == copy
        res.cpus == copy.cpus
        res.memory == copy.memory
        res.disk == copy.disk

    }

    def 'should create resources obj given the cluster config' () {

        given:
        def cfg = new ClusterConfig([maxCpus: 10, maxMemory: "4GB", maxDisk: "100GB" ], 'worker' )

        when:
        def res1 = new Resources(cfg)
        then:
        res1.cpus == 10
        res1.memory == MemoryUnit.of('4GB')
        res1.disk == MemoryUnit.of('100GB')

        when:
        def res2 = new Resources(new ClusterConfig([:], 'master'))
        then:
        res2.cpus == SysHelper.getAvailCpus()
        res2.memory.toGiga() == SysHelper.getAvailMemory().toGiga()
        res2.disk.toGiga() == SysHelper.getAvailDisk().toGiga()

    }

    def 'should validate task holder object' () {

        given:
        def holder = new TaskHolder( task: Mock(IgBaseTask), worker: UUID.randomUUID(), submitTimestamp: System.currentTimeMillis(), started: true, startTimestamp: System.currentTimeMillis()+100)
        when:
        def copy = holder.clone()
        then:
        copy == holder
        copy.task == holder.task
        copy.worker == holder.worker
        copy.submitTimestamp == holder.submitTimestamp
        copy.started == holder.started
        copy.startTimestamp == holder.startTimestamp

        when:
        copy.completed = true
        then:
        true

    }

    def 'should create node data object' () {

        given:
        def config = new ClusterConfig([:])
        when:
        def data = Protocol.NodeData.create(config)
        then:
        data.hostName == SysHelper.getHostName()
        data.resources.cpus == Runtime.runtime.availableProcessors()
        data.resources.memory > MemoryUnit.of(0)
        data.resources.disk > MemoryUnit.of(0)
        data.bootTimeMillis > 0

    }

    def 'should validate node idle method' () {

        given:
        def config = new ClusterConfig([:])
        def data = Protocol.NodeData.create(config)
        data.idleTimestamp = ts

        expect:
        data.isIdle(to) == expected

        where:
        expected    | to                    | ts
        false       | Duration.of('1 m')    | System.currentTimeMillis() - 100
        true        | Duration.of('1 m')    | System.currentTimeMillis() - 100_000

    }


}
