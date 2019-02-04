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
        cmd.printDefault0(config, buffer)

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
        cmd.printProperties0(config, buffer)

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

        cmd.printCanonical0(config, buffer)
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
        cmd.printFlatten0(config, buffer)

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

        cmd.printFlatten0(config, buffer)
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

    def 'should handle variables' () {
        given:
        def folder = Files.createTempDirectory('test')
        def CONFIG = folder.resolve('nextflow.config')

        CONFIG.text = '''
        X1 = 'SOMETHING'
        X2 = [X1]
        X3 = [p:111, q:'bbb']
        
        params {
          alpha = ["${X1}/a", "b", "c"]
          delta = [ X2, 'z' ]
          gamma = [p: "${X1}/a", q: X2]
          omega = X3 
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
        buffer.toString() == '''\
        X1 = 'SOMETHING'
        X2 = ['SOMETHING']
        X3 = [p:111, q:'bbb']
        
        params {
           alpha = ['SOMETHING/a', 'b', 'c']
           delta = [['SOMETHING'], 'z']
           gamma = [p:'SOMETHING/a', q:['SOMETHING']]
           omega = [p:111, q:'bbb']
        }
        '''.stripIndent()

        when:
        def result = new ConfigSlurper().parse(buffer.toString())
        then:
        result.X1 == 'SOMETHING'
        result.X2 == ['SOMETHING']
        result.X3 == [p:111, q:'bbb']
        result.params.alpha == ['SOMETHING/a', 'b', 'c']
        result.params.delta == [['SOMETHING'], 'z']
        result.params.gamma == [p:'SOMETHING/a', q:['SOMETHING']]
        result.params.omega == [p:111, q:'bbb']

    }

}
