/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
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
