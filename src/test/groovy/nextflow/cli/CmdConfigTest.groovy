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

import spock.lang.Specification
import nextflow.cli.CmdConfig.OrderedProperties
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CmdConfigTest extends Specification {

    def 'should default notation' () {

        given:
        def config = new ConfigObject()
        config.process.executor = 'slurm'
        config.process.queue = 'long'
        config.docker.enabled = true

        final buffer = new ByteArrayOutputStream()
        def cmd = new CmdConfig()

        when:
        cmd.printDefault(config, buffer)

        then:
        buffer.toString() == '''
                process {
                \texecutor='slurm'
                \tqueue='long'
                }
                docker.enabled=true
                '''
                .stripIndent().leftTrim()

    }

    def 'should properties notation' () {

        given:
        def mem = { return 0 }
        def config = new ConfigObject()
        config.process.executor = 'slurm'
        config.process.queue = 'long'
        config.process.memory = mem
        config.docker.enabled = true
        config.process.omega = "Hi' there"

        final buffer = new ByteArrayOutputStream()
        def cmd = new CmdConfig()

        when:
        cmd.printProperties(config, buffer)

        then:
        buffer.toString() == """
                docker.enabled=true
                process.executor=slurm
                process.memory=${mem.toString()}
                process.omega=Hi' there
                process.queue=long
                """
                .stripIndent().leftTrim()

    }
    def 'should canonical notation' () {

        given:
        ByteArrayOutputStream buffer
        ConfigObject config
        def cmd = new CmdConfig()

        when:
        buffer = new ByteArrayOutputStream()

        config = new ConfigObject()
        config.process.executor = 'slurm'
        config.process.queue = 'long'
        config.docker.enabled = true
        config.dummy = new ConfigObject() // <-- empty config object should not be print
        config.mail.from = 'yo@mail.com'
        config.mail.smtp.host = 'mail.com'
        config.mail.smtp.port = 25
        config.mail.smtp.user = 'yo'

        cmd.printCanonical(config, buffer)
        then:
        buffer.toString() == '''
                    docker {
                       enabled = true
                    }

                    mail {
                       from = 'yo@mail.com'
                       smtp {
                          host = 'mail.com'
                          port = 25
                          user = 'yo'
                       }
                    }

                    process {
                       executor = 'slurm'
                       queue = 'long'
                    }
                    '''
                    .stripIndent().leftTrim()

    }

    def 'should flatten notation' () {

        given:
        ByteArrayOutputStream buffer
        ConfigObject config
        def cmd = new CmdConfig()

        when:
        buffer = new ByteArrayOutputStream()
        config = new ConfigObject()
        config.process.executor = 'slurm'
        config.process.queue = 'long'
        config.docker.enabled = true
        cmd.printFlatten(config, buffer)

        then:
        buffer.toString() == '''
                docker.enabled = true
                process.executor = 'slurm'
                process.queue = 'long'
                '''
                .stripIndent().leftTrim()

        when:
        buffer = new ByteArrayOutputStream()
        config = new ConfigObject()
        config.foo = "Hi' there"

        cmd.printFlatten(config, buffer)
        then:
        buffer.toString() == "foo = 'Hi\\' there'\n"


    }


    def 'should sort property keys' () {

        given:
        def props = new OrderedProperties()
        props.setProperty('omega', '1')
        props.setProperty('beta',  '3')
        props.setProperty('delta', '2')
        props.setProperty('alpha', '4')


        when:
        def e = props.keys()
        then:
        e.nextElement() == 'alpha'
        e.nextElement() == 'beta'
        e.nextElement() == 'delta'
        e.nextElement() == 'omega'

    }


    def 'should create from a properties object' () {

        given:
        def config = new Properties()
        config.'omega' = 1
        config.'alpha' = 4
        config.'delta.y' = 'Hello'
        config.'delta.z' = 'world'

        when:
        def props = new OrderedProperties(config)
        then:
        props.'omega' == 1
        props.'alpha' == 4
        props.'delta.y' == 'Hello'
        props.'delta.z' == 'world'
        props.keys().toSet() == ['alpha', 'delta.y', 'delta.z','omega'] as Set
    }

}
