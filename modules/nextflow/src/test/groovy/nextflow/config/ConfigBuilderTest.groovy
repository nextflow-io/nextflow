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

package nextflow.config

import java.nio.file.Files
import java.nio.file.Path

import nextflow.SysEnv
import nextflow.exception.AbortOperationException
import nextflow.exception.ConfigParseException
import nextflow.extension.FilesEx
import nextflow.secret.SecretsLoader
import nextflow.util.ConfigHelper
import spock.lang.Ignore
import spock.lang.Specification
import spock.lang.Unroll
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ConfigBuilderTest extends Specification {

    def 'build config object' () {

        setup:
        def env = [PATH:'/local/bin', HOME:'/home/my']
        def builder = [:] as ConfigBuilder

        when:
        def config = builder.build(env,null)

        then:
        ('PATH' in config.env )
        ('HOME' in config.env )
        !('XXX' in config.env )

        config.env.PATH == '/local/bin'
        config.env.HOME == '/home/my'

    }

    def 'build config object 2' () {

        setup:
        def builder = [:] as ConfigBuilder
        def env = [HOME:'/home/my', PATH:'/local/bin', 'dot.key.name':'any text']

        def text1 = '''
        task { field1 = 1; field2 = 'hola'; }
        env { alpha = 'a1'; beta  = 'b1'; HOME="$HOME:/some/path"; }
        '''

        def text2 = '''
        task { field2 = 'Hello' }
        env { beta = 'b2'; delta = 'd2'; HOME="$HOME:/local/path"; XXX="$PATH:/xxx"; YYY = "$XXX:/yyy"; WWW = "${WWW?:''}:value"   }
        '''

        when:
        def config1 = builder.build(env, [text1])
        def config2 = builder.build(env, [text1, text2])

        // note: configuration object can be modified like any map
        config2.env ['ZZZ'] = '99'

        then:
        config1.task.field1 == 1
        config1.task.field2 == 'hola'
        config1.env.HOME == '/home/my:/some/path'
        config1.env.PATH == '/local/bin'
        config1.env.'dot.key.name' == 'any text'

        config2.task.field1 == 1
        config2.task.field2 == 'Hello'
        config2.env.alpha == 'a1'
        config2.env.beta == 'b2'
        config2.env.delta == 'd2'
        config2.env.HOME == '/home/my:/local/path'
        config2.env.XXX == '/local/bin:/xxx'
        config2.env.PATH == '/local/bin'
        config2.env.YYY == '/local/bin:/xxx:/yyy'
        config2.env.ZZZ == '99'
        config2.env.WWW == ':value'

    }

    def 'build config object 3' () {

        setup:
        def builder = [:] as ConfigBuilder
        def env = [HOME:'/home/my', PATH:'/local/bin', 'dot.key.name':'any text']

        def text1 = '''
        task { field1 = 1; field2 = 'hola'; }
        env { alpha = 'a1'; beta  = 'b1'; HOME="$HOME:/some/path"; }
        params { demo = 1  }
        params.test = 2
        '''

        when:
        def config1 = builder.build(env, [text1])

        then:
        config1.task.field1 == 1
        config1.task.field2 == 'hola'
        config1.env.HOME == '/home/my:/some/path'
        config1.env.PATH == '/local/bin'
        config1.env.'dot.key.name' == 'any text'

        config1.params.test == 2
        config1.params.demo == 1

    }

    def 'build config object 4' () {

        setup:
        def builder = [:] as ConfigBuilder
        builder.baseDir = Path.of('/base/path')

        def text = '''
        params.p = "$baseDir/1"
        params {
            q = "$baseDir/2"
            x = "$projectDir/3"
            y = "$launchDir/4"
            z = "$outputDir/5"
        }
        '''

        when:
        def cfg = builder.build([text])
        then:
        cfg.params.p == '/base/path/1'
        cfg.params.q == '/base/path/2'
        cfg.params.x == '/base/path/3'
        cfg.params.y == "${Path.of('.').toRealPath()}/4"
        cfg.params.z == "${Path.of('results').complete()}/5"

    }

    def 'should check for a valid profile' () {

        given:
        def builder = new ConfigBuilder()

        when:
        builder.profile = 'omega'
        builder.checkValidProfile([])
        then:
        def ex = thrown(AbortOperationException)
        ex.message == "Unknown configuration profile: 'omega'"

        when:
        builder.profile = 'omigo'
        builder.checkValidProfile(['alpha','delta','omega'])
        then:
        ex = thrown(AbortOperationException)
        ex.message == """
            Unknown configuration profile: 'omigo'

            Did you mean one of these?
                omega

            """
            .stripIndent().leftTrim()

        when:
        builder.profile = 'delta'
        builder.checkValidProfile(['alpha','delta','omega'])
        then:
        true

        when:
        builder.profile = 'delta,omega'
        builder.checkValidProfile(['alpha','delta','omega'])
        then:
        true

        when:
        builder.profile = 'delta,bravo'
        builder.checkValidProfile(['alpha','delta','omega'])
        then:
        ex = thrown(AbortOperationException)
        ex.message == "Unknown configuration profile: 'bravo'"

    }

    def 'should warn about missing attribute' () {

        given:
        def file = Files.createTempFile('test','config')
        file.deleteOnExit()
        file.text =
                '''
                params.foo = HOME
                '''

        when:
        def cfg = new ConfigBuilder().build([file])
        then:
        cfg.params.foo == System.getenv('HOME')

        when:
        file.text =
                '''
                params.foo = bar
                '''
        new ConfigBuilder().build([file])
        then:
        def e = thrown(ConfigParseException)
        e.message == "Unknown config attribute `bar` -- check config file: ${file.toRealPath()}".toString()

    }

    def 'should render missing variables' () {
        given:
        def file = Files.createTempFile('test',null)

        file.text =
                '''
                foo = 'xyz'
                bar = "$SCRATCH/singularity_images_nextflow"
                '''

        when:
        def builder = new ConfigBuilder()
                .setShowMissingVariables(true)
        def cfg = builder.build([file])
        def str = ConfigHelper.toCanonicalString(cfg)
        then:
        str == '''\
            foo = 'xyz'
            bar = '$SCRATCH/singularity_images_nextflow'
            '''.stripIndent()

        and:
        builder.warnings[0].startsWith('Unknown config attribute `SCRATCH`')
        cleanup:
        file?.delete()
    }

    def 'should report fully qualified missing attribute'  () {

        given:
        def file = Files.createTempFile('test','config')
        file.deleteOnExit()

        when:
        file.text =
                '''
                params.x = foo.bar
                '''
        new ConfigBuilder().build([file])
        then:
        def e = thrown(ConfigParseException)
        e.message == "Unknown config attribute `foo.bar` -- check config file: ${file.toRealPath()}".toString()

    }

    def 'should collect config files' () {

        given:
        def slurper = ConfigParserFactory.create()
        def file1 = Files.createTempFile('test1', null)
        def file2 = Files.createTempFile('test2', null)
        def result = new ConfigObject()
        def builder = new ConfigBuilder()

        file1.text = 'foo = 1'
        file2.text = 'bar = 2'

        when:
        builder.merge0(result, slurper, file1)
        builder.merge0(result, slurper, file2)
        then:
        result.foo == 1
        result.bar == 2
        builder.parsedConfigFiles == [file1, file2]

        cleanup:
        file1?.delete()
        file2?.delete()
    }

    def 'should merge profiles' () {
        given:
        def ENV = [:]
        def result
        def builder = new ConfigBuilder()

        def CONFIG = '''
                process.container = 'base'
                process.executor = 'local'
                
                profiles {
                    cfg1 {
                      process.executor = 'sge'
                      process.queue = 'short'
                    }

                    cfg2 {
                        process.executor = 'batch'
                        process.queue = 'long'
                    }

                    docker {
                        docker.enabled = true
                        process.container = 'foo/1'
                    }

                    singularity {
                        singularity.enabled = true
                        process.queue = 'cn-el7'
                        process.container = 'bar-2.img'
                    }
                }
                '''

        when:
        builder.profile = 'cfg1'
        result = builder.build(ENV, [CONFIG])
        then:
        result.process.container == 'base'
        result.process.executor == 'sge'
        result.process.queue == 'short'

        when:
        builder.profile = 'cfg2'
        result = builder.build(ENV, [CONFIG])
        then:
        result.process.container == 'base'
        result.process.executor == 'batch'
        result.process.queue == 'long'

        when:
        builder.profile = 'cfg1,docker'
        result = builder.build(ENV, [CONFIG])
        then:
        result.process.container ==  'foo/1'
        result.process.executor == 'sge'
        result.process.queue == 'short'
        result.docker.enabled
        !result.isSet('singularity')

        when:
        builder.profile = 'cfg1,singularity'
        result = builder.build(ENV, [CONFIG])
        then:
        result.process.container ==  'bar-2.img'
        result.process.executor == 'sge'
        result.process.queue == 'cn-el7'
        result.singularity.enabled
        !result.isSet('docker')

        when:
        builder.profile = 'cfg2,singularity'
        result = builder.build(ENV, [CONFIG])
        then:
        result.process.container ==  'bar-2.img'
        result.process.executor == 'batch'
        result.process.queue == 'cn-el7'
        result.singularity.enabled
        !result.isSet('docker')

        when:
        builder.profile = 'missing'
        builder.build(ENV, [CONFIG])
        then:
        thrown(AbortOperationException)
    }

    def 'should return all profiles' () {
        given:
        def ENV = [:]
        def result
        def builder = new ConfigBuilder()

        def CONFIG = '''

                profiles {
                    cfg1 {
                      process.executor = 'sge'
                      process.queue = 'short'
                    }

                    cfg2 {
                        process.executor = 'batch'
                        process.queue = 'long'
                    }

                    docker {
                        docker.enabled = true
                        process.container = 'foo/1'
                    }

                    singularity {
                        singularity.enabled = true
                        process.queue = 'cn-el7'
                        process.container = 'bar-2.img'
                    }
                }
                '''

        when:
        builder.showAllProfiles = true
        result = builder.build(ENV, [CONFIG])
        then:
        result.profiles.cfg1.process.executor == 'sge'
        result.profiles.cfg1.process.queue == 'short'
        result.profiles.cfg2.process.executor == 'batch'
        result.profiles.cfg2.process.queue == 'long'
        result.profiles.docker.docker.enabled
        result.profiles.docker.process.container == 'foo/1'
        result.profiles.singularity.singularity.enabled
        result.profiles.singularity.process.queue == 'cn-el7'
        result.profiles.singularity.process.container == 'bar-2.img'
    }

    def 'should resolve process withLabel inside a profile' () {

        given:
        def folder = Files.createTempDirectory('test')
        def file1 = folder.resolve('test.conf')
        file1.text = '''
            process {
                cpus = 2

                withName: bar {
                    cpus = 4
                }

                withLabel: foo {
                    cpus = 8
                }
            }
            '''

        when:
        def cfg1 = new ConfigBuilder().build([file1])
        then:
        cfg1.process.cpus == 2
        cfg1.process.'withName:bar'.cpus == 4
        cfg1.process.'withLabel:foo'.cpus == 8
        cfg1.process == [cpus: 2, 'withName:bar': [cpus:4], 'withLabel:foo': [cpus: 8]]

        // make sure the same config is returned when the file is included
        when:
        def file2 = folder.resolve('nextflow.config')
        file2.text = """
            includeConfig "$file1" 
            """

        def cfg2 = new ConfigBuilder().build([file2])
        then:
        cfg1 == cfg2

        cleanup:
        folder?.deleteDir()

    }

    def 'should resolve ext config' () {

        given:
        def folder = Files.createTempDirectory('test')
        def file1 = folder.resolve('test.conf')
        file1.text = '''
            process {
                ext { args = "Hello World!" } 
                cpus = 1 
                withName:BAR {
                    ext { args = "Ciao mondo!" } 
                    cpus = 2
                }
            }
            '''

        when:
        def cfg1 = new ConfigBuilder().build([file1])
        then:
        cfg1.process.cpus == 1
        cfg1.process.ext.args == 'Hello World!'
        cfg1.process.'withName:BAR'.cpus == 2
        cfg1.process.'withName:BAR'.ext.args == "Ciao mondo!"

        cleanup:
        folder?.deleteDir()
    }

    // issue 2422 - https://github.com/nextflow-io/nextflow/issues/2422
    // ideally this should behave as the previous test
    @Ignore
    def 'should resolve ext config with properties' () {

        given:
        def folder = Files.createTempDirectory('test')
        def file1 = folder.resolve('test.conf')
        file1.text = '''
            process {
                ext.args = "Hello World!" 
                cpus = 1 
                withName:BAR {
                    ext.args = "Ciao mondo!"
                    cpus = 2
                }
            }
            '''

        when:
        def cfg1 = new ConfigBuilder().build([file1])
        then:
        cfg1.process.cpus == 1
        cfg1.process.ext.args == 'Hello World!'
        cfg1.process.'withName:BAR'.cpus == 2
        cfg1.process.'withName:BAR'.ext.args == "Ciao mondo!"

        cleanup:
        folder?.deleteDir()
    }

    def 'should access top params from profile' () {
        given:
        def folder = Files.createTempDirectory('test')
        def file1 = folder.resolve('file1.conf')

        file1.text = """
            params.alpha = 1 
            params.delta = 2 
            
            profiles {
                foo {
                    params.delta = 20 
                    params.gamma = 30 
                    
                    process {  
                        cpus = params.alpha
                    }
                }
            }
            """

        when:
        def cfg = new ConfigBuilder().setProfile('foo').build([file1])
        then:
        cfg.process == [cpus: 1]
        cfg.params == [alpha: 1, delta: 20, gamma: 30]

        cleanup:
        folder?.deleteDir()
    }

    def 'should access top params from profile [2]' () {
        given:
        def folder = Files.createTempDirectory('test')
        def file1 = folder.resolve('file1.conf')

        file1.text = """
            params {
                alpha = 1 
                delta = 2 
            }
            
            profiles {
                foo {
                    params {
                        delta = 20 
                        gamma = 30
                    }

                    process {  
                        cpus = params.alpha
                    }
                }
            }
            """

        when:
        def cfg = new ConfigBuilder().setProfile('foo').build([file1])
        then:
        cfg.process == [cpus: 1]
        cfg.params == [alpha: 1, delta: 20, gamma: 30]

        cleanup:
        folder?.deleteDir()
    }

    def 'should merge params two profiles' () {
        given:
        def folder = Files.createTempDirectory('test')
        def file1 = folder.resolve('file1.conf')

        file1.text = '''    
            profiles {
                foo {
                    params {
                        alpha = 1 
                        delta = 2 
                    }
                }

                bar {
                    params {
                        delta = 20 
                        gamma = 30 
                    }
                    
                    process {
                        cpus = params.alpha
                    }                
                }
            }
            '''

        when:
        def cfg = new ConfigBuilder().setProfile('foo,bar') .build([file1])
        then:
        cfg.process.cpus == 1
        cfg.params.alpha == 1
        cfg.params.delta == 20
        cfg.params.gamma == 30

        cleanup:
        folder?.deleteDir()
    }

    def 'should merge profiles with conditions' () {
        given:
        def folder = Files.createTempDirectory("mergeprofiles")
        def main = folder.resolve('main.conf')
        def test = folder.resolve('test.conf')
        def process = folder.resolve('process.conf')

        main.text = '''
        params {
            load_config = null
            present = true
        }
        
        profiles {
            test { includeConfig 'test.conf' }
        }
        
        if (params.load_config) {
            includeConfig 'process.conf'
        }        
        '''
        test.text = '''
        params {
            load_config = true
        }    
        '''

        process.text = '''
        process {
            withName: FOO {
                ext.args = '--quiet'
            }
        }        
        params{
            another = true
        }
        '''

        when:
        def cfg = new ConfigBuilder().setProfile('test').build([main])
        then:
        cfg.process.'withName:FOO'
        cfg.params.load_config == true
        cfg.params.present == true
        cfg.params.another == true

        cleanup:
        folder?.deleteDir()
    }

    def 'should build config object with secrets' () {
        given:
        SecretsLoader.instance.reset()
        and:
        def folder = Files.createTempDirectory('test')
        def secrets  = folder.resolve('store.json')
        and:
        secrets.text = '''
            [
              {
                "name": "FOO",
                "value": "ciao"
              }
            ]
            '''
        and:
        FilesEx.setPermissions(secrets, 'rw-------')
        SysEnv.push(NXF_SECRETS_FILE:secrets.toAbsolutePath().toString())
        and:
        def text = '''
        params.p = "$baseDir/1"
        params.s = "$secrets.FOO/2"
        '''

        when:
        def cfg = new ConfigBuilder().setBaseDir(folder).build([text])
        then:
        cfg.params.p == "$folder/1"
        cfg.params.s == 'ciao/2'

        cleanup:
        folder?.deleteDir()
        SysEnv.pop()
    }

    def 'should build include config with secrets' () {
        given:
        SecretsLoader.instance.reset()
        and:
        def folder = Files.createTempDirectory('test')
        def secrets  = folder.resolve('store.json')
        and:
        secrets.text = '''
            [
              {
                "name": "ALPHA",
                "value": "one"
              },
              {
                "name": "DELTA",
                "value": "two"
              },
              {
                "name": "GAMMA",
                "value": "three"
              }
            ]
            '''
        FilesEx.setPermissions(secrets, 'rw-------')
        SysEnv.push(NXF_SECRETS_FILE:secrets.toAbsolutePath().toString())
        and:

        def configMain = folder.resolve('nextflow.config')
        def snippet1 = folder.resolve('config1.txt')
        def snippet2 = folder.resolve('config2.txt')
        and:

        configMain.text = """
        p1 = secrets.ALPHA
        includeConfig "$snippet1"
        """

        snippet1.text = """
        p2 = secrets.DELTA
        includeConfig "$snippet2" 
        """

        snippet2.text = '''
        p3 = secrets.GAMMA
        p4 = 'some value'
        '''

        when:
        def config = new ConfigBuilder().setBaseDir(folder).build([configMain])
        then:
        config.p1 == 'one'
        config.p2 == 'two'
        config.p3 == 'three'
        config.p4 == 'some value'

        cleanup:
        folder?.deleteDir()
        SysEnv.pop()
    }

    def 'should build config with secrets in the include path' () {
        given:
        SecretsLoader.instance.reset()
        and:
        def folder = Files.createTempDirectory('test')
        def configMain = folder.resolve('nextflow.config')
        def snippet1 = folder.resolve('config1.txt')
        def snippet2 = folder.resolve('config2.txt')

        and:
        def secrets  = folder.resolve('store.json')
        and:
        secrets.text = """
            [
              {
                "name": "SECRET_FILE1",
                "value": "${snippet1.toAbsolutePath()}"
              },
              {
                "name": "SECRET_FILE2",
                "value": "${snippet2.toAbsolutePath()}"
              }
            ]
            """
        FilesEx.setPermissions(secrets, 'rw-------')
        SysEnv.push(NXF_SECRETS_FILE:secrets.toAbsolutePath().toString())
        and:
        configMain.text = '''
        p1 = 'one'
        // include config via a secret property
        includeConfig secrets.SECRET_FILE1
        '''

        snippet1.text = '''
        p2 = 'two'
        // include config via string interpolation using a secret 
        includeConfig "$secrets.SECRET_FILE2"
        '''

        snippet2.text = '''
        p3 = 'three'
        '''

        when:
        def config = new ConfigBuilder().setBaseDir(folder).build([configMain])
        then:
        config.p1 == 'one'
        config.p2 == 'two'
        config.p3 == 'three'

        cleanup:
        folder?.deleteDir()
        SysEnv.pop()
    }

    def 'should not render secret values' () {
        given:
        SecretsLoader.instance.reset()
        and:
        def folder = Files.createTempDirectory('test')
        def configMain = folder.resolve('nextflow.config')
        def snippet1 = folder.resolve('config1.txt')
        def snippet2 = folder.resolve('config2.txt')

        and:
        def secrets  = folder.resolve('store.json')
        and:
        secrets.text = """
            [
              {
                "name": "SECRET_FILE1",
                "value": "${snippet1.toAbsolutePath()}"
              },
              {
                "name": "SECRET_FILE2",
                "value": "${snippet2.toAbsolutePath()}"
              },
              {
                "name": "ALPHA",
                "value": "one"
              },
              {
                "name": "DELTA",
                "value": "two"
              },
              {
                "name": "GAMMA",
                "value": "three"
              }
            ]
            """
        FilesEx.setPermissions(secrets, 'rw-------')
        SysEnv.push(NXF_SECRETS_FILE:secrets.toAbsolutePath().toString())
        and:
        configMain.text = '''
        p1 = secrets.ALPHA
        foo.p1 = secrets.ALPHA 
        bar {
          p1 = secrets.ALPHA 
        }
        // include config via a secret property
        includeConfig secrets.SECRET_FILE1
        '''

        snippet1.text = '''
        p2 = "$secrets.DELTA"
        foo.p2 = "$secrets.DELTA"
        bar {
          p2 = "$secrets.DELTA"
        }
        // include config via string interpolation using a secret 
        includeConfig "$secrets.SECRET_FILE2"
        '''

        snippet2.text = '''
        p3 = "$secrets.GAMMA"
        foo.p3 = "$secrets.GAMMA"
        bar {
          p3 = "$secrets.GAMMA"
        }
        '''

        when:
        def config = new ConfigBuilder()
                .setBaseDir(folder)
                .setStripSecrets(true)
                .build([configMain])
        then:
        config.p1 == 'secrets.ALPHA'
        config.p2 == 'secrets.DELTA'
        config.p3 == 'secrets.GAMMA'
        and:
        config.foo.p1 == 'secrets.ALPHA'
        config.foo.p2 == 'secrets.DELTA'
        config.foo.p3 == 'secrets.GAMMA'
        and:
        config.bar.p1 == 'secrets.ALPHA'
        config.bar.p2 == 'secrets.DELTA'
        config.bar.p3 == 'secrets.GAMMA'

        cleanup:
        folder?.deleteDir()
        SysEnv.pop()
    }

}
