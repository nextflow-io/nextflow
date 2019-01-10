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

package nextflow.cloud

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class LauncherConfigTest extends Specification{

    def 'should render config object' () {

        given:
        def config = [
                imageId: 'ami-43f49030',
                instanceType: 't2.micro',
                subnetId: 'subnet-05222a43',
                driver: 'aws',
                keyHash: 'ssh-rsa AAxxxx',
                userName: 'pditommaso',
                autoscale: [
                        enabled: true,
                        maxInstances: 5,
                        terminateWhenIdle: true,
                ]
        ]


        def cloud = new CloudConfig(config) .setClusterName('my-cluster')
        def scale = cloud.getAutoscale()

        when:
        def str1 = cloud.toConfigObject().prettyPrint()
        then:
        str1 == '''
                imageId='ami-43f49030'
                instanceType='t2.micro'
                subnetId='subnet-05222a43'
                driver='aws'
                keyHash='ssh-rsa AAxxxx'
                userName='pditommaso'
                autoscale {
                \tenabled=true
                \tmaxInstances=5
                \tterminateWhenIdle=true
                }
                clusterName='my-cluster'
                '''
                .stripIndent().leftTrim()

        when:
        def str2 = cloud.toCloudConfigObject().prettyPrint()
        then:
        str2 == '''
                cloud {
                \timageId='ami-43f49030'
                \tinstanceType='t2.micro'
                \tsubnetId='subnet-05222a43'
                \tdriver='aws'
                \tkeyHash='ssh-rsa AAxxxx'
                \tuserName='pditommaso'
                \tautoscale {
                \t\tenabled=true
                \t\tmaxInstances=5
                \t\tterminateWhenIdle=true
                \t}
                \tclusterName='my-cluster'
                }
                '''
                .stripIndent().leftTrim()

        when:
        def str3 = scale.toConfigObject().prettyPrint()
        then:
        str3 == '''
                enabled=true
                maxInstances=5
                terminateWhenIdle=true
                '''
                .stripIndent().leftTrim()

        when:
        def str4 = scale.renderCloudConfigObject()
        then:
        str4 == cloud.renderCloudConfigObject()
        str4 == str2
    }

}
