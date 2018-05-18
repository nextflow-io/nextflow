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

import spock.lang.Specification
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
        cmd.sort = true
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
                    process {
                       executor = 'slurm'
                       queue = 'long'
                    }
                    
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
                process.executor = 'slurm'
                process.queue = 'long'
                docker.enabled = true
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


    def 'should parse config file' () {
        given:
        def folder = Files.createTempDirectory('test')
        def CONFIG = folder.resolve('nextflow.config')

        CONFIG.text = '''
        manifest {
            author = 'me'
            mainScript = 'foo.nf'
        }
        
        process {
          queue = 'cn-el7'
          cpus = 4 
          memory = 10.GB
          time = 5.h
          ext.other = { 10.GB * task.attempt }
        }
        '''
        def buffer = new ByteArrayOutputStream()
        // command definition 
        def cmd = new CmdConfig()
        cmd.launcher = new Launcher(options: new CliOptions(config: [CONFIG.toString()]))
        cmd.stdout = buffer
        cmd.args = [ '.' ]

        when:
        cmd.run()
        
        then:
        buffer.toString() == '''
        manifest {
           author = 'me'
           mainScript = 'foo.nf'
        }
        
        process {
           queue = 'cn-el7'
           cpus = 4
           memory = '10 GB'
           time = '5h'
           ext {
              other = { 10.GB * task.attempt }
           }
        }
        '''
        .stripIndent().leftTrim()

        cleanup:
        folder.deleteDir()
    }


}
