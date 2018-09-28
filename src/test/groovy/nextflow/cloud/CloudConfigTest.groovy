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

package nextflow.cloud
import java.nio.file.Files

import nextflow.Const
import nextflow.config.ConfigBuilder
import nextflow.util.Duration
import nextflow.util.MemoryUnit
import spock.lang.Specification
import spock.util.environment.RestoreSystemProperties
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CloudConfigTest extends Specification {


    def 'should populate config object' () {

        given:
        def config = [
                imageId: 'ami-123',
                instanceType: 'r3.abc',
                securityGroup: 'sg-df72b9ba',
                keyName: 'eurokey',
                subnetId: 'subnet-05222a43',
                bootStorageSize: '10GB',
                instanceStorageDevice: '/dev/abc',
                instanceStorageMount: '/scratch',
                sharedStorageId: 'fs-1803efd1',
                sharedStorageMount: '/mnt/efs',
                spotPrice: 0.55,
                instanceRole: 'foo-role',
                dockerPull: 'cbcrg/image:tag'
        ]

        when:
        def cloud = new CloudConfig(config)
        then:
        with(cloud) {
            imageId == 'ami-123'
            instanceType == 'r3.abc'
            securityGroups == ['sg-df72b9ba']
            keyName == 'eurokey'
            bootStorageSize == MemoryUnit.of('10 GB')
            instanceStorageDevice == '/dev/abc'
            instanceStorageMount == '/scratch'
            sharedStorageId == 'fs-1803efd1'
            sharedStorageMount == '/mnt/efs'
            spotPrice == '0.55'
            dockerPull == ['cbcrg/image:tag']
            instanceRole == 'foo-role'
        }

    }


    def 'should return a pretty formatted string' () {

        given:
        def config = [
                imageId: 'ami-123',
                instanceType: 'r3.abc',
                securityGroup: 'sg-df72b9ba',
                keyHash: 'ssh-rsa AAAAB3Nza...',
                subnetId: 'subnet-05222a43',
                bootStorageSize: '10GB',
                spotPrice: 0.55
        ]

        when:
        def text = new CloudConfig(config).prettyPrint()
        then:
        text == """
                - bootStorageSize: '10GB'
                - imageId: 'ami-123'
                - instanceType: 'r3.abc'
                - keyHash: 'ssh-rsa AAAAB3Nza...'
                - securityGroup: 'sg-df72b9ba'
                - spotPrice: 0.55
                - subnetId: 'subnet-05222a43'
                """
                .stripIndent().leftTrim()


    }

    def 'should read the config file and render the config' () {

        given:
        def folder = Files.createTempDirectory('test')
        def cfgFile = folder.resolve('nxf.config')
        cfgFile.text = '''

            cloud {
              imageId = 'ami-43f49030'
              instanceType = 't2.micro'
              subnetId = 'subnet-05222a43'
              sharedStorageId = 'fs-1803efd1'
              sharedStorageMount = '/mnt/efs'

              autoscale {
                  enabled = true
                  instanceType = 'm3.xlarge'
                  starvingTimeout = '6 min'
                  terminateWhenIdle = true
                  spotPrice = 0.35
              }
            }
        '''

        when:
        def config = new ConfigBuilder().buildGivenFiles(cfgFile)
        def cloud = CloudConfig.create(config)

        then:
        cloud.imageId == 'ami-43f49030'
        cloud.instanceType == 't2.micro'
        cloud.subnetId == 'subnet-05222a43'
        cloud.sharedStorageId == 'fs-1803efd1'
        cloud.sharedStorageMount ==  '/mnt/efs'
        with(cloud.autoscale) {
            enabled == true
            instanceType == 'm3.xlarge'
            starvingTimeout == Duration.of('6 min')
            terminateWhenIdle == true
            spotPrice == '0.35'
        }

        expect:
        cloud.prettyPrint() ==  '''
                                - imageId: 'ami-43f49030'
                                - instanceType: 't2.micro'
                                - sharedStorageId: 'fs-1803efd1'
                                - sharedStorageMount: '/mnt/efs'
                                - subnetId: 'subnet-05222a43'
                                - autoscale:
                                  - enabled: true
                                  - instanceType: 'm3.xlarge'
                                  - spotPrice: 0.35
                                  - starvingTimeout: '6 min'
                                  - terminateWhenIdle: true
                                '''
                                .stripIndent().leftTrim()
        cleanup:
        folder?.deleteDir()

    }



    def 'should set instanceType and imageId attributes' () {

        given:
        def config = [
                imageId: 'ami-123',
                instanceType: 'r3.abc',
                keyName: 'eurokey',
        ]

        when:
        def cloud1 = new CloudConfig(config)
                                .setRole('master')
                                .setImageId('ami-89fs89fd')
                                .setInstanceType('m2.xlarge')
                                .build()
        then:
        cloud1.imageId == 'ami-89fs89fd'
        cloud1.instanceType == 'm2.xlarge'
        cloud1.keyName == 'eurokey'
        cloud1.role == 'master'

    }


    def 'should convert to a config object' () {

        given:
        def config = [
                imageId: 'ami-123',
                instanceType: 'r3.abc',
                securityGroup: 'sg-df72b9ba',
                keyName: 'eurokey',
                subnetId: 'subnet-05222a43',
                bootStorageSize: '10GB',
                spotPrice: 0.55
        ]

        expect:
        new CloudConfig(config).renderCloudConfigObject() ==
                '''
                cloud {
                	imageId='ami-123'
                	instanceType='r3.abc'
                	securityGroup='sg-df72b9ba'
                	keyName='eurokey'
                	subnetId='subnet-05222a43'
                	bootStorageSize='10GB'
                	spotPrice=0.55
                }
                '''
                .stripIndent().leftTrim()
    }


    def 'should convert to a config object with scopes' () {

        given:
        def config =  [
                imageId: 'ami-123',
                instanceType: 'r3.abc',
        ]

        when:
        def configObj1 = new CloudConfig(config) .renderCloudConfigObject()
        then:
        configObj1 ==  '''
                        cloud {
                        	imageId='ami-123'
                        	instanceType='r3.abc'
                        }
                        '''
                        .stripIndent().leftTrim()



    }

    def 'should set property keyName and userName' () {
        given:
        def config = [
                imageId: 'ami-123',
                instanceType: 'r3.abc',
                keyName: 'cloud-provided-key'
        ]

        //
        // when is specified the attribute `keyName`
        // the user is not created in the instance
        // and the keyFile and keyHash are null
        //
        when:
        def cloud = new CloudConfig(config) .build()
        then:
        cloud.keyName == 'cloud-provided-key'
        cloud.keyFile == null
        cloud.privateKeyFile == null
        cloud.keyHash == null
        !cloud.createUser

    }

    def 'should set keyFile and userName' () {
        given:
        def key = Files.createTempFile('nf', 'test')
        key.text = 'ssh-public-key'

        and:
        def config = [
                imageId: 'ami-123',
                keyFile: key,
                userName: 'foo'
        ]

        when:
        def cloud = new CloudConfig(config) .build()
        then:
        cloud.userName == 'foo'
        cloud.keyFile == key
        cloud.keyHash == 'ssh-public-key'

        cleanup:
        key.delete()

    }

    @RestoreSystemProperties
    def 'should set keyFile and local userName' () {

        given:
        final HOME = Files.createTempDirectory('test')
        System.setProperty('user.home', HOME.toString())
        HOME.resolve('.ssh').mkdir()
        HOME.resolve('.ssh/id_rsa.pub').text = 'ssh-rsa AAAAB3NzaC1yc2EAAAABIwAAAQEA0fPOx7mcUAl3YIhqlllMOZsZMpXKjzfCNYK'
        HOME.resolve('.ssh/id_rsa').text = '--- private key here ---'

        def config = [
                imageId: 'ami-123',
                instanceType: 'r3.abc',
                autoscale: [
                        spotPrice: 0.33
                ]
        ]

        //
        // when the keyName is not specified it reads the local key file
        // and creates a user with the same name as the local user name
        //
        when:
        def cloud = new CloudConfig(config)
                            .setClusterName('my-cloud')
                            .setRole('master')
                            .build()
        then:
        cloud.keyName == null
        cloud.keyFile == HOME.resolve('.ssh/id_rsa.pub')
        cloud.privateKeyFile == HOME.resolve('.ssh/id_rsa')
        cloud.keyHash == 'ssh-rsa AAAAB3NzaC1yc2EAAAABIwAAAQEA0fPOx7mcUAl3YIhqlllMOZsZMpXKjzfCNYK'
        cloud.createUser
        cloud.role == 'master'
        cloud.clusterName == 'my-cloud'

        when:
        def copy = cloud.getAutoscale()
        then:
        copy.keyName == null
        copy.keyHash == 'ssh-rsa AAAAB3NzaC1yc2EAAAABIwAAAQEA0fPOx7mcUAl3YIhqlllMOZsZMpXKjzfCNYK'
        copy.createUser
        copy.clusterName == 'my-cloud'
        copy.role == 'worker' // <-- alwasy dy definition

        cleanup:
        HOME?.deleteDir()

    }


    def 'should parse autoscale data' () {

        when:
        def cloud = new CloudConfig([
                imageId: 'ami-43f49030',
                instanceType: 't2.micro',
                spotPrice: '0.05'
        ])
        then:
        !cloud.autoscale.enabled

        when:
        cloud = new CloudConfig([
                imageId: 'ami-43f49030',
                instanceType: 't2.micro',
                spotPrice: '0.05',
                autoscale: true
        ])
        then:
        cloud.instanceType == 't2.micro'
        with(cloud.autoscale) {
            enabled

            instanceType == 't2.micro'
            imageId == 'ami-43f49030'
            spotPrice == '0.05'
        }

        when:
        cloud = new CloudConfig([
                imageId: 'ami-43f49030',
                instanceType: 't2.micro',
                spotPrice: '0.05',
                autoscale: [
                        imageId: 'ami-89ds89',
                        instanceType: 'm3.large',
                        spotPrice: '0.33',
                        starvingTimeout: '2 m',
                        terminateWhenIdle: true
                ]
        ])
        then:
        cloud.instanceType == 't2.micro'
        with(cloud.autoscale) {
            enabled
            instanceType == 'm3.large'
            imageId == 'ami-89ds89'
            spotPrice == '0.33'
            terminateWhenIdle
            starvingTimeout == Duration.of('2 min')
        }
    }

    def 'should create a autoscale config object from a map' () {

        when:
        def auto = CloudConfig.Autoscale.create([:])
        then:
        auto.enabled
        !auto.terminateWhenIdle
        auto.starvingTimeout == Duration.of('5 min')

        when:
        auto = CloudConfig.Autoscale.create([
                enabled:true,
                imageId: 'ami-89ds',
                instanceType:'i-32424354',
                spotPrice: '0.50',
                terminateWhenIdle: true,
                starvingTimeout: '1 m',
        ])
        then:
        auto.enabled
        auto.terminateWhenIdle
        auto.imageId == 'ami-89ds'
        auto.instanceType == 'i-32424354'
        auto.spotPrice == '0.50'
        auto.starvingTimeout == Duration.of('1 min')

        when:
        CloudConfig.Autoscale.create([unknown: true])
        then:
        thrown(IllegalArgumentException)
    }

    def 'should return nextflow config obj' () {

        when:
        def cloud1 = new CloudConfig([
                imageId: 'ami-43f49030',
                instanceType: 't2.micro',
        ])
        then:
        cloud1.nextflow.version == Const.APP_VER
        cloud1.nextflow.mode == 'ignite'

        when:
        def cloud2 = new CloudConfig([
                imageId: 'ami-43f49030',
                instanceType: 't2.micro',
                nextflow: [
                        version: '1.2.0',
                        pull: 'pditommaso/hello',
                        options: '-Dx=y',
                        trace: 'com.foo.bar'
                ]
        ])

        then:
        cloud2.nextflow.version == '1.2.0'
        cloud2.nextflow.pull == 'pditommaso/hello'
        cloud2.nextflow.options == '-Dx=y'
        cloud2.nextflow.mode == 'ignite'
        cloud2.nextflow.trace == 'com.foo.bar'

    }

    def 'should throw an exception when specifying unknown fields' ()  {

        given:
        def cfg = [enabled: true, starvingTimeout: '5min']

        when:
        def scale = new CloudConfig.Autoscale(cfg, null)
        then:
        scale.enabled
        scale.starvingTimeout == Duration.of('5min')

        when:
        new CloudConfig.Autoscale([enabled: true, unknown: '5min'], null)
        then:
        thrown(IllegalArgumentException)

    }

    def 'should create an autoscale config object' () {

        when:
        def scale1 = CloudConfig.Autoscale.create(null)
        then:
        !scale1.enabled

        when:
        def scale2 = CloudConfig.Autoscale.create(true)
        then:
        with( scale2 ) {
            enabled
            minInstances == 1
            maxInstances == Integer.MAX_VALUE
            role == 'worker'
        }

        when:
        def scale3 = CloudConfig.Autoscale.create([imageId: 'abc', minInstances: 10, maxInstances: 100])
        then:
        with(scale3) {
            enabled
            imageId == 'abc'
            minInstances == 10
            maxInstances == 100
            role == 'worker'
        }

    }
}
