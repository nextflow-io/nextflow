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

package nextflow.cloud.azure.config

import com.azure.compute.batch.models.BatchNodeCommunicationMode
import com.azure.compute.batch.models.ImageVerificationType
import com.azure.compute.batch.models.OSType
import com.google.common.hash.Hashing
import nextflow.util.CacheHelper
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
        !opts.virtualMachineImageId
        opts.osType == OSType.LINUX
        opts.verification == ImageVerificationType.VERIFIED
        opts.targetCommunicationMode == BatchNodeCommunicationMode.SIMPLIFIED
    }

    def 'should configure a custom compute gallery image' () {
        when:
        def opts = new AzPoolOpts([
            virtualMachineImageId: '/subscriptions/abc/resourceGroups/rg/providers/Microsoft.Compute/galleries/g/images/d/versions/1.0.0',
            sku: 'batch.node.ubuntu 22.04',
            verification: 'unverified',
            osType: 'linux',
        ])
        then:
        opts.virtualMachineImageId == '/subscriptions/abc/resourceGroups/rg/providers/Microsoft.Compute/galleries/g/images/d/versions/1.0.0'
        opts.sku == 'batch.node.ubuntu 22.04'
        opts.verification == ImageVerificationType.UNVERIFIED
        opts.osType == OSType.LINUX
    }

    def 'should parse the verification value' () {
        expect:
        new AzPoolOpts([verification: VALUE]).verification == EXPECTED
        where:
        VALUE        | EXPECTED
        null         | ImageVerificationType.VERIFIED
        'verified'   | ImageVerificationType.VERIFIED
        'unverified' | ImageVerificationType.UNVERIFIED
        'any'        | null
    }

    def 'should reject an invalid verification value' () {
        when:
        new AzPoolOpts([verification: 'bogus'])
        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains('verification')
    }

    def 'should parse the osType value' () {
        expect:
        new AzPoolOpts([osType: VALUE]).osType == EXPECTED
        where:
        VALUE     | EXPECTED
        null      | OSType.LINUX
        'linux'   | OSType.LINUX
        'windows' | OSType.WINDOWS
    }

    def 'should reject an invalid osType value' () {
        when:
        new AzPoolOpts([osType: 'macos'])
        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains('osType')
    }

    def 'should parse the target communication mode' () {
        expect:
        new AzPoolOpts([targetCommunicationMode: VALUE]).targetCommunicationMode == EXPECTED
        where:
        VALUE        | EXPECTED
        null         | BatchNodeCommunicationMode.SIMPLIFIED
        'simplified' | BatchNodeCommunicationMode.SIMPLIFIED
        'classic'    | BatchNodeCommunicationMode.CLASSIC
        'SIMPLIFIED' | BatchNodeCommunicationMode.SIMPLIFIED
        'CLASSIC'    | BatchNodeCommunicationMode.CLASSIC
    }

    def 'should reject an invalid target communication mode' () {
        when:
        new AzPoolOpts([targetCommunicationMode: 'bogus'])
        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains('targetCommunicationMode')
    }

    private static String hash(AzPoolOpts opts) {
        opts.funnel(Hashing.murmur3_128().newHasher(), CacheHelper.HashMode.STANDARD).hash().toString()
    }

    def 'pool hash should differ when image config differs' () {
        given:
        def base = new AzPoolOpts()
        def gallery = new AzPoolOpts([virtualMachineImageId: '/subscriptions/x/resourceGroups/rg/providers/Microsoft.Compute/galleries/g/images/d/versions/1'])
        def unverified = new AzPoolOpts([verification: 'unverified'])
        expect:
        hash(base) != hash(gallery)
        hash(base) != hash(unverified)
    }

    def 'pool hash should differ when communication mode differs' () {
        given:
        def base = new AzPoolOpts()
        def classic = new AzPoolOpts([targetCommunicationMode: 'classic'])
        expect:
        hash(base) != hash(classic)
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
