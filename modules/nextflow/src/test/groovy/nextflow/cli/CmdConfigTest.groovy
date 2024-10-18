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
 */

package nextflow.cli

import java.nio.file.Files

import nextflow.SysEnv
import nextflow.exception.AbortOperationException
import nextflow.extension.FilesEx
import nextflow.plugin.Plugins
import nextflow.secret.SecretsLoader
import spock.lang.IgnoreIf
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CmdConfigTest extends Specification {

    def cleanup() {
        Plugins.stop()
    }

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

    def 'should print the value of a config option' () {

        given:
        def cmd = new CmdConfig()
        and:
        def config = new ConfigObject()
        config.process.executor = 'slurm'
        config.process.queue = 'long'
        config.docker.enabled = true
        and:
        def buffer

        when:
        buffer = new ByteArrayOutputStream()
        cmd.printValue0(config, 'process.executor', buffer)
        then:
        buffer.toString() == 'slurm\n'

        when:
        buffer = new ByteArrayOutputStream()
        cmd.printValue0(config, 'does.not.exist', buffer)
        then:
        def e = thrown(AbortOperationException)
        e.message == "Configuration option 'does.not.exist' not found"

    }

    def 'should print config using json' () {
        given:
        ByteArrayOutputStream buffer
        ConfigObject config
        def cmd = new CmdConfig(sort: true)

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

        cmd.printJson0(config, buffer)
        then:
        buffer.toString() == '''\
                    {
                        "docker": {
                            "enabled": true
                        },
                        "mail": {
                            "from": "yo@mail.com",
                            "smtp": {
                                "host": "mail.com",
                                "port": 25,
                                "user": "yo"
                            }
                        },
                        "process": {
                            "executor": "slurm",
                            "queue": "long"
                        }
                    }
                    '''
                    .stripIndent()
    }

    def 'should print config using yaml' () {
        given:
        ByteArrayOutputStream buffer
        ConfigObject config
        def cmd = new CmdConfig(sort: true)

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

        cmd.printYaml0(config, buffer)
        then:
        buffer.toString() == '''\
                    docker:
                      enabled: true
                    mail:
                      from: yo@mail.com
                      smtp:
                        host: mail.com
                        port: 25
                        user: yo
                    process:
                      executor: slurm
                      queue: long
                    '''
                    .stripIndent()
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
          gamma = [p: "${X1}/a", q: X2, 'r-r': 'X1']
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
           gamma = [p:'SOMETHING/a', q:['SOMETHING'], 'r-r':'X1']
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
        result.params.gamma == [p:'SOMETHING/a', q:['SOMETHING'], 'r-r': 'X1']
        result.params.omega == [p:111, q:'bbb']

    }


    @IgnoreIf({System.getenv('NXF_SMOKE')})
    def 'should resolve remote config' () {
        given:
        def buffer = new ByteArrayOutputStream()
        def cmd = new CmdConfig(
                args: ['https://github.com/nextflow-io/hello'],
                showAllProfiles: true,
                launcher: Mock(Launcher),
                stdout: buffer  )

        when:
        cmd.run()
        then:
        buffer.toString() == '''\
            process {
               container = 'quay.io/nextflow/bash'
            }
            '''
            .stripIndent()
    }

    @IgnoreIf({System.getenv('NXF_SMOKE')})
    def 'should resolve profiles into profiles config' () {
        given:
        def folder = Files.createTempDirectory('test')
        def CONFIG = folder.resolve('nextflow.config')

        CONFIG.text = '''
            params {
               foo = 'baz'
            }
            
            profiles {
               test {
                  params {
                    foo = 'foo'
                  }
                  profiles {                      
                      debug {
                        cleanup = false
                      }                    
                  }
               }
            }        
        '''

        def buffer = new ByteArrayOutputStream()
        // command definition
        def cmd = new CmdConfig(showAllProfiles: true)
        cmd.launcher = new Launcher(options: new CliOptions(config: [CONFIG.toString()]))
        cmd.stdout = buffer
        cmd.args = ['.']

        when:
        cmd.run()
        def result = new ConfigSlurper().parse(buffer.toString())

        then:
        result.params.foo == 'baz'
        result.profiles.test.params.foo == 'foo'
        result.profiles.debug.cleanup == false
    }


    def 'should remove secrets for config' () {
        given:
        SecretsLoader.instance.reset()
        and:
        def folder = Files.createTempDirectory('test')
        and:
        def secrets  = folder.resolve('store.json')
        and:
        secrets.text = "[ ]"
        FilesEx.setPermissions(secrets, 'rw-------')
        SysEnv.push(NXF_SECRETS_FILE:secrets.toAbsolutePath().toString())

        and:
        def CONFIG = folder.resolve('nextflow.config')
        CONFIG.text = '''
        process { 
            queue = secrets.MYSTERY
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
        process {
           queue = 'secrets.MYSTERY'
        }
        '''
            .stripIndent().leftTrim()

        cleanup:
        folder?.deleteDir()
        SysEnv.pop()
    }

    def 'should render json config output' () {
        given:
        def folder = Files.createTempDirectory('test')
        def CONFIG = folder.resolve('nextflow.config')

        CONFIG.text = '''
        manifest {
            author = 'me'
            mainScript = 'foo.nf'
        }
        
        process {
          cpus = 4 
          queue = 'cn-el7'
          memory = { 10.GB }   
          ext.other = { 10.GB * task.attempt }
        }
        '''
        def buffer = new ByteArrayOutputStream()
        // command definition
        def cmd = new CmdConfig(outputFormat: 'json')
        cmd.launcher = new Launcher(options: new CliOptions(config: [CONFIG.toString()]))
        cmd.stdout = buffer
        cmd.args = [ '.' ]

        when:
        cmd.run()

        then:
        buffer.toString() == '''
        {
            "manifest": {
                "author": "me",
                "mainScript": "foo.nf"
            },
            "process": {
                "cpus": 4,
                "queue": "cn-el7",
                "memory": "{ 10.GB }",
                "ext": {
                    "other": "{ 10.GB * task.attempt }"
                }
            }
        }
        '''
            .stripIndent().leftTrim()

        cleanup:
        folder.deleteDir()
    }

    def 'should render yaml config output' () {
        given:
        def folder = Files.createTempDirectory('test')
        def CONFIG = folder.resolve('nextflow.config')

        CONFIG.text = '''
        manifest {
            author = 'me'
            mainScript = 'foo.nf'
        }
        
        process {
          cpus = 4 
          queue = 'cn-el7'
          memory = { 10.GB }   
          ext.other = { 10.GB * task.attempt }
        }
        '''
        def buffer = new ByteArrayOutputStream()
        // command definition
        def cmd = new CmdConfig(outputFormat: 'yaml')
        cmd.launcher = new Launcher(options: new CliOptions(config: [CONFIG.toString()]))
        cmd.stdout = buffer
        cmd.args = [ '.' ]

        when:
        cmd.run()

        then:
        buffer.toString() == '''
            manifest:
              author: me
              mainScript: foo.nf
            process:
              cpus: 4
              queue: cn-el7
              memory: '{ 10.GB }'
              ext:
                other: '{ 10.GB * task.attempt }'
            '''
                .stripIndent().leftTrim()

        cleanup:
        folder.deleteDir()
    }
}
