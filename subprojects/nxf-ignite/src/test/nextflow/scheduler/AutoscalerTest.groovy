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

import nextflow.cloud.CloudConfig.Autoscale
import nextflow.cloud.CloudDriver
import nextflow.cloud.types.CloudInstanceStatus
import nextflow.cloud.types.CloudInstanceType
import nextflow.executor.IgBaseTask
import nextflow.processor.TaskId
import nextflow.scheduler.Protocol.NodeData
import nextflow.scheduler.Protocol.Resources
import nextflow.scheduler.Protocol.TaskHolder
import nextflow.scheduler.Protocol.TaskResources
import nextflow.util.Duration
import nextflow.util.MemoryUnit
import spock.lang.Specification
import spock.lang.Unroll
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AutoscalerTest extends Specification {

    def 'should limit the number of killed nodes' () {

        given:
        def scaler = Spy(Autoscaler)
        scaler.getClusterSize() >> 10
        def list = [ Mock(NodeData),Mock(NodeData),Mock(NodeData),Mock(NodeData),Mock(NodeData) ]

        when:
        scaler.scalerConfig = Autoscale.create([:])
        def killList1 = scaler.ensureMinSize( list )
        then:
        killList1.size() == 5

        when:
        scaler.scalerConfig = Autoscale.create([minInstances: 8])
        def killList2 = scaler.ensureMinSize( list )
        then:
        killList2.size() == 2
    }


    def 'should invoke aws client api' () {

        given:
        def driver = Mock(CloudDriver)
        def config = Autoscale.create([instanceType: 'm4.xlarge'])
        def instances = ['i1','i2','i3','i4','i5']
        def scaler = new Autoscaler(driver: driver, scalerConfig: config)
        scaler.enabled = true

        when:
        scaler.requestNewInstances( 5 )
        then:
        !scaler.enabled
        scaler.pendingInstanceIds as List == instances

        then:
        1 * driver.launchInstances(5, config) >> instances

        then:
        1 * driver.waitInstanceStatus(instances, CloudInstanceStatus.STARTED)

        then:
        driver.tagInstances(instances, config)

    }

    def 'should infer the number of instances to launch' () {

        given:
        def config = Autoscale.create([instanceType: 'm4.xlarge'])
        def driver = Mock(CloudDriver)
        def scaler = Spy(Autoscaler)
        scaler.driver = driver
        scaler.scalerConfig = config
        scaler.workerNodes = [:]
        def instanceType = new CloudInstanceType(id:'m4.xlarge', cpus: 4, memory: MemoryUnit.of('16GB'))

        // requesting: 50 cpus
        // instance type: 4 cpus
        // ==> 50 / 4 = 12.5 ==> 13 instances requested
        when:
        scaler.requestNewCpus(50)
        then:
        1 * driver.describeInstanceType('m4.xlarge') >> instanceType
        1 * scaler.requestNewInstances(13)
        1 * driver.launchInstances(13, config) >> []

        // requesting: 50 cpus
        // instance type: 4 cpus
        // ==> 13 instances needed
        // max cluster size: 10 instances
        // ==> final request: 10 instances
        when:
        scaler.scalerConfig = Autoscale.create([instanceType: 'm4.xlarge', maxInstances: 10])
        scaler.requestNewCpus(50)
        then:
        1 * driver.describeInstanceType('m4.xlarge') >> instanceType
        1 * scaler.requestNewInstances(10)
        1 * driver.launchInstances(10, _) >> []

        // requesting: 50 cpus
        // instance type: 4 cpus
        // ==> 13 cpus needed
        // max cluster size: 10
        // current cluster size: 3 instances
        // ==> final request: 7 instances
        when:
        scaler.scalerConfig = Autoscale.create([instanceType: 'm4.xlarge', maxInstances: 10])
        scaler.requestNewCpus(50)
        then:
        1 * scaler.getClusterSize() >> 3
        1 * driver.describeInstanceType('m4.xlarge') >> instanceType
        1 * scaler.requestNewInstances(7)
        1 * driver.launchInstances(7, _) >> []

    }

    def 'should resize the cluster if needed'() {

        given:
        def config = Autoscale.create([instanceType: 'm4.xlarge', minInstances: 10])
        def driver = Mock(CloudDriver)
        def scaler = Spy(Autoscaler)
        scaler.driver = driver
        scaler.scalerConfig = config
        scaler.enabled = true

        when:
        scaler.checkClusterSize()
        then:
        1 * scaler.getClusterSize() >> 6
        1 * scaler.requestNewInstances(4)
        1 * driver.launchInstances(4,_) >> []

    }

    @Unroll
    def 'should validate instance type resources #test' () {
        given:
        def config = Autoscale.create([instanceType: 'm4.xlarge', minInstances: 10])
        def driver = Mock(CloudDriver)
        def scaler = Spy(Autoscaler)
        scaler.driver = driver
        scaler.scalerConfig = config
        scaler.enabled = true

        def task = Mock(IgBaseTask)
        task.getTaskId() >> { TaskId.of(100) }
        task.getResources() >> { new Protocol.TaskResources(cpus: req_cpus, memory: req_mem ? MemoryUnit.of(req_mem) : null) }

        when:
        def result = scaler.canRunOnNewInstance(task, 'm4.xlarge')
        then:
        1 * driver.describeInstanceType('m4.xlarge') >> { new CloudInstanceType(cpus: tot_cpus, memory: MemoryUnit.of(tot_mem))  }
        result == expected

        where:
        test    | req_cpus     | req_mem   | tot_cpus  | tot_mem    | expected
        1       | 4            | null      | 4         | '8GB'      | true
        2       | 5            | null      | 4         | '8GB'      | false
        3       | 1            | '8GB'     | 2         | '8GB'      | true
        4       | 1            | '9GB'     | 2         | '8GB'      | false

    }


    def "should check if there's a not fulfilling the resource request of a given task" () {

        given:
        def _1_GB  =  MemoryUnit.of('1GB')
        def _4_GB  =  MemoryUnit.of('4GB')
        def _8_GB  =  MemoryUnit.of('8GB')
        def _16_GB  =  MemoryUnit.of('16GB')

        def topology = [
                'x1': new NodeData(free: new Resources(cpus: 0, memory: _4_GB), resources: new Resources(cpus: 4, memory:_4_GB)),
                'x2': new NodeData(free: new Resources(cpus: 8, memory: _4_GB), resources: new Resources(cpus: 8, memory:_8_GB)),
                "x3": new NodeData(free: new Resources(cpus: 6, memory: _8_GB), resources: new Resources(cpus: 18, memory:_16_GB))
        ]

        def config = Autoscale.create([instanceType: 'm4.xlarge', minInstances: 10])
        def driver = Mock(CloudDriver)
        def scaler = Spy(Autoscaler)
        scaler.driver = driver
        scaler.scalerConfig = config
        scaler.enabled = true
        scaler.workerNodes = topology


        def task = Mock(IgBaseTask)
        task.getTaskId() >> { TaskId.of(taskId) }
        task.getResources() >> { new Protocol.TaskResources(cpus: req_cpus, memory: req_mem ? MemoryUnit.of(req_mem) : null) }

        expect:
        scaler.canRunOnExistingNodes(task) == expected

        where:
        taskId  | req_cpus  | req_mem   | expected
        1       | 1         | null      | 1     // can execute
        2       | 8         | null      | 1     // can execute
        3       | 10        | null      | 0     // cant
        4       | 8         | '8GB'     | 0     // cant
        5       | 6         | '8GB'     | 1     // can execute
        6       | 2         | '16GB'    | 0     // cant
        7       | 32        | '2GB'     | 2     // resources overflow

    }

    private createTaskMock( Map params ) {

        def now = System.currentTimeMillis()
        def task = Mock(IgBaseTask);
        task.getTaskId() >> TaskId.of(params.id)
        task.getResources() >> new TaskResources(cpus: params.cpus as int, memory: params.memory as MemoryUnit)

        def holder = new TaskHolder(task)

        if( params.start ) {
            holder.started = true
            holder.startTimestamp = now
        }

        return holder
    }

    @Unroll
    def 'should request new instances for starving task #task_id' () {

        given:
        def _4_GB  =  MemoryUnit.of('4GB')
        def _8_GB  =  MemoryUnit.of('8GB')
        def _16_GB  =  MemoryUnit.of('16GB')

        def topology = [
                'x2': new NodeData(free: new Resources(cpus: 8, memory: _4_GB), resources: new Resources(cpus: 8, memory:_16_GB)),
                "x3": new NodeData(free: new Resources(cpus: 6, memory: _8_GB), resources: new Resources(cpus: 8, memory:_16_GB))
        ]

        def config = Autoscale.create([instanceType: 'm4.xlarge', minInstances: 10, starvingTimeout: timeout as Duration])
        def driver = Mock(CloudDriver)
        def scaler = Spy(Autoscaler)
        scaler.driver = driver
        scaler.scalerConfig = config
        scaler.workerNodes = topology
        scaler.enabled = true

        def holder = createTaskMock(id: task_id, cpus: cpus, memory: memory, start: start)

        when:
        scaler.scheduledTasks = [(holder.task.taskId): holder]
        scaler.waitingTimestamp = timestamp ? (System.currentTimeMillis() - Duration.of(timestamp).millis) : 0
        scaler.checkStarvingTasks()
        then:
        _ * driver.describeInstanceType('m4.xlarge') >> new CloudInstanceType(cpus: 16, memory: _16_GB)
        make_request * scaler.requestNewCpus(cpus)
        make_request * driver.launchInstances(1,_) >> []

        where:
        task_id | cpus  | memory | start   | timeout | timestamp | make_request
        1       | 8     | '8 GB' | false   | '5 min' | '6 min'   | 1             // not enough free resources ==> request new instance
        2       | 8     | '8 GB' | false   | '5 min' | '1 min'   | 0             // not enough free resources BUT too recent (timestamp < timeout) ==> NO request new instance
        3       | 8     | '8 GB' | true    | '5 min' | '6 min'   | 0             // not enough free resources BUT already started ==> NO request new instance
        4       | 16    | '16GB' | false   | '5 min' | '6 min'   | 1             // no node in the topology to launch this node ==> launch a new instance
        5       | 100   | '16GB' | false   | '5 min' | '6 min'   | 0             // no instance available ==> no launch instances for this task

    }


    def 'validate check idle nodes' () {

        given:
        def now = System.currentTimeMillis()
        def _10_mins_ago = now - Duration.of('10 min').toMillis()
        def _57_mins_ago = now - Duration.of('57 min').toMillis()
        def _90_mins_ago = now - Duration.of('90 min').toMillis()
        def n1 = UUID.fromString('11111111-9b9f-4b45-bebc-aae51fad9d18')
        def n2 = UUID.fromString('22222222-2b3b-4a98-b771-f20d1ed50c83')
        def n3 = UUID.fromString('33333333-663c-4760-8300-4d0ff90ef7cd')
        def n4 = UUID.fromString('44444444-663c-4760-8300-4d0ff90ef7cd')
        def n5 = UUID.fromString('55555555-663c-4760-8300-4d0ff90ef7cd')

        def topology = [
                (n1): new NodeData(nodeId: n1, instanceId: 'i1', bootTimeMillis: _10_mins_ago, idleTimestamp: now-10),
                (n2): new NodeData(nodeId: n2, instanceId: 'i2', bootTimeMillis: _10_mins_ago, idleTimestamp: now-10_000),
                (n3): new NodeData(nodeId: n3, instanceId: 'i3', bootTimeMillis: _10_mins_ago, idleTimestamp: now-10_000),
                (n4): new NodeData(nodeId: n4, instanceId: 'i4', bootTimeMillis: _57_mins_ago, idleTimestamp: now-10_000),
                (n5): new NodeData(nodeId: n5, instanceId: 'i5', bootTimeMillis: _90_mins_ago, idleTimestamp: now-10_000),
        ]

        def config = Autoscale.create([terminateWhenIdle: true, idleTimeout: '1 s', terminationPolicy: 'eager'])
        def driver = Mock(CloudDriver)
        def scaler = Spy(Autoscaler)
        scaler.driver = driver
        scaler.scalerConfig = config
        scaler.workerNodes = topology
        scaler.enabled = true

        when:
        scaler.checkIdleNodes()
        then:
        1 * scaler.getLocalNodeId() >> n3
        1 * scaler.applyTerminationPolicy([n2,n4,n5])
        1 * scaler.killNodes( ['i2','i4','i5'] )
        1 * scaler.ensureMinSize( ['i2','i4','i5'] )
        1 * scaler.requestTerminateInstances(['i2','i4','i5'])
        1 * driver.terminateInstances(['i2','i4','i5'])
        1 * driver.waitInstanceStatus(['i2','i4','i5'], CloudInstanceStatus.TERMINATED)

        // scaler not enabled
        // ==> no nodes killed
        when:
        scaler.enabled = false
        scaler.checkIdleNodes()
        then:
        0 * scaler.getLocalNodeId()
        0 * scaler.applyTerminationPolicy(_)
        0 * scaler.killNodes(_)

        // `terminateWhenIdle` config attribute is not set
        // ==> not nodes killed
        when:
        scaler.enabled = true
        scaler.scalerConfig = Autoscale.create([terminateWhenIdle: false, idleTimeout: '1 s'])
        scaler.checkIdleNodes()
        then:
        0 * scaler.getLocalNodeId()
        0 * scaler.applyTerminationPolicy(_)
        0 * scaler.killNodes(_)

//        // termination policy is not `eager`
//        // ==> kill only node up between 55 to 60 minutes (i.e. node # 4)
//        when:
//        scaler.enabled = true
//        scaler.scalerConfig = Autoscale.create([terminateWhenIdle: true, idleTimeout: '1 s' ])
//        scaler.checkIdleNodes()
//        then:
//        1 * scaler.getLocalNodeId() >> n3
//        1 * scaler.applyTerminationPolicy([n2,n4,n5])
//        1 * scaler.killNodes(['i4'])

    }

}
