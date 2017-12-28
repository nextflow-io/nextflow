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
