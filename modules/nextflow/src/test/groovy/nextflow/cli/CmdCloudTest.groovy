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

package nextflow.cli

import java.nio.file.Files

import nextflow.cloud.CloudDriver
import nextflow.exception.AbortOperationException
import spock.lang.Specification
import spock.util.environment.RestoreSystemProperties

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CmdCloudTest extends Specification {

    def 'should validate the cluster name' () {

        given:
        def cmd = new CmdCloud()

        when:
        cmd.checkName('hello')
        then:
        true

        when:
        cmd.checkName('Hello-1_2')
        then:
        true

        when:
        cmd.checkName('Hello world')
        then:
        thrown(AbortOperationException)

        when:
        cmd.checkName(null)
        then:
        thrown(AbortOperationException)

        when:
        cmd.checkName('')
        then:
        thrown(AbortOperationException)

        when:
        cmd.checkName('123')
        then:
        thrown(AbortOperationException)

        when:
        cmd.checkName('1abc')
        then:
        thrown(AbortOperationException)

    }

    @RestoreSystemProperties
    def 'should create a cloud config object' () {
        given:
        def folder = Files.createTempDirectory('test')
        folder.resolve('.ssh').mkdir()
        folder.resolve('.ssh/id_rsa.pub').text = 'ssh-rsa fake'
        System.setProperty('user.home', folder.toString())

        def cmd = Spy(CmdCloud)
        cmd.config = [ cloud: [
            imageId: 'ami-blah',
            instanceType: 'm4.xxlarge'
        ]]
        cmd.driver = Mock(CloudDriver)

        when:
        def config1 = cmd.makeConfig('my-cloud')
        then:
        with(config1) {
            clusterName == 'my-cloud'
            imageId == 'ami-blah'
            instanceType == 'm4.xxlarge'
            spotPrice == null
            driverName == null
        }

        when:
        cmd.imageId = 'ami-09f9gd'
        cmd.instanceType = 'c3.small'
        cmd.spotPrice = '0.43'
        cmd.driverName = 'gcloud'
        def config2 = cmd.makeConfig('your-cloud')
        then:
        with(config2) {
            clusterName == 'your-cloud'
            imageId == 'ami-09f9gd'
            instanceType == 'c3.small'
            spotPrice == '0.43'
            driverName == 'gcloud'
        }

        cleanup:
        folder?.deleteDir()
    }



}
