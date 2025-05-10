/*
 * Copyright 2013-2024, Seqera Labs
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
 *
 */

package nextflow.cloud.azure.config

import nextflow.util.Duration
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AzPoolOptsTest extends Specification {

    def 'should create pool options' () {
        when:
        def opts = new AzPoolOpts()
        then:
        !opts.runAs
        !opts.privileged
        opts.publisher == AzPoolOpts.DEFAULT_PUBLISHER
        opts.offer == AzPoolOpts.DEFAULT_OFFER
        opts.sku == AzPoolOpts.DEFAULT_SKU
        opts.vmType == AzPoolOpts.DEFAULT_VM_TYPE
        opts.fileShareRootPath == '/mnt/batch/tasks/fsmounts'
        opts.vmCount == 1
        !opts.autoScale
        !opts.scaleFormula
        !opts.schedulePolicy
        opts.scaleInterval == AzPoolOpts.DEFAULT_SCALE_INTERVAL
        opts.maxVmCount == opts.vmCount *3
        !opts.registry
        !opts.userName
        !opts.password
        !opts.virtualNetwork
        !opts.lowPriority
        !opts.startTask.script
        !opts.startTask.privileged
    }

    def 'should create pool with custom options' () {
        when:
        def opts = new AzPoolOpts([
            runAs:'foo',
            privileged: true,
            publisher: 'some-pub',
            offer: 'some-offer',
            sku: 'some-sku',
            vmType: 'some-vmtype',
            vmCount: 10,
            autoScale: true,
            scaleFormula: 'some-formula',
            schedulePolicy: 'some-policy',
            scaleInterval: Duration.of('10s'),
            maxVmCount: 100,
            registry: 'some-reg',
            userName: 'some-user',
            password: 'some-pwd',
            virtualNetwork: 'some-vnet',
            lowPriority: true,
            startTask: [
                script: 'echo hello-world',
                privileged: true
            ]
        ])
        then:
        opts.runAs == 'foo'
        opts.privileged
        opts.publisher == 'some-pub'
        opts.offer == 'some-offer'
        opts.sku == 'some-sku'
        opts.vmType == 'some-vmtype'
        opts.fileShareRootPath == ''
        opts.vmCount == 10
        opts.autoScale
        opts.scaleFormula == 'some-formula'
        opts.schedulePolicy == 'some-policy'
        opts.scaleInterval == Duration.of('10s')
        opts.maxVmCount == 100
        opts.registry == 'some-reg'
        opts.userName == 'some-user'
        opts.password == 'some-pwd'
        opts.virtualNetwork == 'some-vnet'
        opts.lowPriority
        opts.startTask.script == 'echo hello-world'
        opts.startTask.privileged
    }

}
