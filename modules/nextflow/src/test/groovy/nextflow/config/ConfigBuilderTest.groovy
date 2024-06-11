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
import java.nio.file.Paths

import nextflow.SysEnv
import nextflow.cli.CliOptions
import nextflow.cli.CmdConfig
import nextflow.cli.CmdNode
import nextflow.cli.CmdRun
import nextflow.cli.Launcher
import nextflow.exception.AbortOperationException
import nextflow.exception.ConfigParseException
import nextflow.extension.FilesEx
import nextflow.secret.SecretsLoader
import nextflow.trace.TraceHelper
import nextflow.util.ConfigHelper
import spock.lang.Ignore
import spock.lang.Specification
import spock.lang.Unroll
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ConfigBuilderTest extends Specification {

    def setup() {
        TraceHelper.testTimestampFmt = '20221001'
    }

    def 'build config object' () {

        setup:
        def env = [PATH:'/local/bin', HOME:'/home/my']
        def builder = [:] as ConfigBuilder

        when:
        def config = builder.buildConfig0(env,null)

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
        def config1 = builder.buildConfig0(env, [text1])
        def config2 = builder.buildConfig0(env, [text1, text2])

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
        def config1 = builder.buildConfig0(env, [text1])

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
        builder.baseDir = Paths.get('/base/path')

        def text = '''
        params.p = "$baseDir/1"
        params {
            q = "$baseDir/2"
            x = "$projectDir/3"
            y = "$launchDir/4"
        }
        '''

        when:
        def cfg = builder.buildConfig0([:], [text])
        then:
        cfg.params.p == '/base/path/1'
        cfg.params.q == '/base/path/2'
        cfg.params.x == '/base/path/3'
        cfg.params.y == "${Paths.get('.').toRealPath()}/4"

    }

    def 'CLI params should override the ones defined in the config file' () {
        setup:
        def file = Files.createTempFile('test',null)
        file.text = '''
        params {
          alpha = 'x'
        }
        params.beta = 'y'
        params.delta = 'Foo'
        params.gamma = params.alpha
        params {
            omega = 'Bar'
        }

        process {
          publishDir = [path: params.alpha]
        }
        '''
        when:
        def opt = new CliOptions()
        def run = new CmdRun(params: [alpha: 'Hello', beta: 'World', omega: 'Last'])
        def result = new ConfigBuilder().setOptions(opt).setCmdRun(run).buildGivenFiles(file)

        then:
        result.params.alpha == 'Hello'  // <-- params defined as CLI options override the ones in the config file
        result.params.beta == 'World'   // <--   as above
        result.params.gamma == 'Hello'  // <--   as above
        result.params.omega == 'Last'
        result.params.delta == 'Foo'
        result.process.publishDir == [path: 'Hello']

        cleanup:
        file?.delete()
    }

    def 'CLI params should override the ones defined in the config file [2]' () {
        setup:
        def file = Files.createTempFile('test',null)
        file.text = '''
        params {
          alpha = 'x'
          beta = 'y'
          delta = 'Foo'
          gamma = params.alpha
          omega = 'Bar'
        }

        process {
          publishDir = [path: params.alpha]
        }
        '''
        when:
        def opt = new CliOptions()
        def run = new CmdRun(params: [alpha: 'Hello', beta: 'World', omega: 'Last'])
        def result = new ConfigBuilder().setOptions(opt).setCmdRun(run).buildGivenFiles(file)

        then:
        result.params.alpha == 'Hello'  // <-- params defined as CLI options override the ones in the config file
        result.params.beta == 'World'   // <--   as above
        result.params.gamma == 'Hello'  // <--   as above
        result.params.omega == 'Last'
        result.params.delta == 'Foo'
        result.process.publishDir == [path: 'Hello']

        cleanup:
        file?.delete()
    }


    def 'CLI params should override the ones in one or more config files' () {
        given:
        def folder = File.createTempDir()
        def configMain = new File(folder,'nextflow.config').absoluteFile
        def snippet1 = new File(folder,'config1.txt').absoluteFile
        def snippet2 = new File(folder,'config2.txt').absoluteFile


        configMain.text = """
        process.name = 'alpha'
        params.one = 'a'
        params.xxx = 'x'
        includeConfig "$snippet1"
        """

        snippet1.text = """
        params.two = 'b'
        params.yyy = 'y'

        process.cpus = 4
        process.memory = '8GB'

        includeConfig("$snippet2")
        """

        snippet2.text = '''
        params.three = 'c'
        params.zzz = 'z'

        process { disk = '1TB' }
        process.resources.foo = 1
        process.resources.bar = 2
        '''

        when:
        def opt = new CliOptions()
        def run = new CmdRun(params: [one: '1', two: 'dos', three: 'tres'])
        def config = new ConfigBuilder().setOptions(opt).setCmdRun(run).buildGivenFiles(configMain.toPath())

        then:
        config.params.one == 1
        config.params.two == 'dos'
        config.params.three == 'tres'
        config.process.name == 'alpha'
        config.params.xxx == 'x'
        config.params.yyy == 'y'
        config.params.zzz == 'z'

        config.process.cpus == 4
        config.process.memory == '8GB'
        config.process.disk == '1TB'
        config.process.resources.foo == 1
        config.process.resources.bar == 2

        cleanup:
        folder?.deleteDir()
    }

    def 'should include config with params' () {
        given:
        def folder = File.createTempDir()
        def configMain = new File(folder,'nextflow.config').absoluteFile
        def snippet1 = new File(folder,'igenomes.config').absoluteFile


        configMain.text = '''
        includeConfig 'igenomes.config'
        '''

        snippet1.text = '''
        params {
          genomes {
            'GRCh37' {
              fasta = "${params.igenomes_base}/genome.fa"
              bwa   = "${params.igenomes_base}/BWAIndex/genome.fa"
            }
          }
        }
        '''

        when:
        def opt = new CliOptions()
        def run = new CmdRun(params: [igenomes_base: 'test'])
        def config = new ConfigBuilder().setOptions(opt).setCmdRun(run).buildGivenFiles(configMain.toPath())

        then:
        config.params.genomes.GRCh37 == [fasta:'test/genome.fa', bwa:'test/BWAIndex/genome.fa']

        cleanup:
        folder?.deleteDir()
    }


    def 'should fetch the config path from env var' () {
        given:
        def folder = File.createTempDir()
        def configMain = new File(folder,'my.config').absoluteFile


        configMain.text = """
        process.name = 'alpha'
        params.one = 'a'
        params.two = 'b'
        """

        // relative path to current dir
        when:
        def config = new ConfigBuilder(env: [NXF_CONFIG_FILE: 'my.config']) .setCurrentDir(folder.toPath()) .build()
        then:
        config.params.one == 'a'
        config.params.two == 'b'
        config.process.name == 'alpha'

        // absolute path
        when:
        config = new ConfigBuilder(env: [NXF_CONFIG_FILE: configMain.toString()]) .build()
        then:
        config.params.one == 'a'
        config.params.two == 'b'
        config.process.name == 'alpha'

        // default should not find it
        when:
        config = new ConfigBuilder() .build()
        then:
        config.params == [:]

        cleanup:
        folder?.deleteDir()
    }

    def 'CLI params should overrides the ones in one or more profiles' () {

        setup:
        def file = Files.createTempFile('test',null)
        file.text = '''
        params.alpha = 'a'
        params.beta = 'b'
        params.delta = 'Foo'
        params.gamma = params.alpha

        params {
            genomes {
                'GRCh37' {
                  bed12   = '/data/genes.bed'
                  bismark = '/data/BismarkIndex'
                  bowtie  = '/data/genome'
                }
            }
        }

        profiles {
          first {
            params.alpha = 'Alpha'
            params.omega = 'Omega'
            params.gamma = 'First'
            process.name = 'Bar'
          }

          second {
            params.alpha = 'xxx'
            params.gamma = 'Second'
            process {
                publishDir = [path: params.alpha]
            }
          }

        }

        '''


        when:
        def opt = new CliOptions()
        def run = new CmdRun(params: [alpha: 'AAA', beta: 'BBB'])
        def config = new ConfigBuilder().setOptions(opt).setCmdRun(run).buildGivenFiles(file)
        then:
        config.params.alpha == 'AAA'
        config.params.beta == 'BBB'
        config.params.delta == 'Foo'
        config.params.gamma == 'AAA'
        config.params.genomes.GRCh37.bed12 == '/data/genes.bed'
        config.params.genomes.GRCh37.bismark == '/data/BismarkIndex'
        config.params.genomes.GRCh37.bowtie  == '/data/genome'

        when:
        opt = new CliOptions()
        run = new CmdRun(params: [alpha: 'AAA', beta: 'BBB'], profile: 'first')
        config = new ConfigBuilder().setOptions(opt).setCmdRun(run).buildGivenFiles(file)
        then:
        config.params.alpha == 'AAA'
        config.params.beta == 'BBB'
        config.params.delta == 'Foo'
        config.params.gamma == 'First'
        config.process.name == 'Bar'
        config.params.genomes.GRCh37.bed12 == '/data/genes.bed'
        config.params.genomes.GRCh37.bismark == '/data/BismarkIndex'
        config.params.genomes.GRCh37.bowtie  == '/data/genome'


        when:
        opt = new CliOptions()
        run = new CmdRun(params: [alpha: 'AAA', beta: 'BBB', genomes: 'xxx'], profile: 'second')
        config = new ConfigBuilder().setOptions(opt).setCmdRun(run).buildGivenFiles(file)
        then:
        config.params.alpha == 'AAA'
        config.params.beta == 'BBB'
        config.params.delta == 'Foo'
        config.params.gamma == 'Second'
        config.params.genomes == 'xxx'
        config.process.publishDir == [path: 'AAA']

        cleanup:
        file?.delete()
    }

    def 'params-file should override params in the config file' () {
        setup:
        def baseDir = Paths.get('/my/base/dir')
        and:
        def params = Files.createTempFile('test', '.yml')
        params.text = '''
            alpha: "Hello" 
            beta: "World" 
            omega: "Last"
            theta: "${baseDir}/something"
            '''.stripIndent()
        and:
        def file = Files.createTempFile('test',null)
        file.text = '''
        params {
          alpha = 'x'
        }
        params.beta = 'y'
        params.delta = 'Foo'
        params.gamma = params.alpha
        params {
            omega = 'Bar'
        }

        process {
          publishDir = [path: params.alpha]
        }
        '''
        when:
        def opt = new CliOptions()
        def run = new CmdRun(paramsFile: params)
        def result = new ConfigBuilder().setOptions(opt).setCmdRun(run).setBaseDir(baseDir).buildGivenFiles(file)

        then:
        result.params.alpha == 'Hello'  // <-- params defined in the params-file overrides the ones in the config file
        result.params.beta == 'World'   // <--   as above
        result.params.gamma == 'Hello'  // <--   as above
        result.params.omega == 'Last'
        result.params.delta == 'Foo'
        result.params.theta == "$baseDir/something"
        result.process.publishDir == [path: 'Hello']

        cleanup:
        file?.delete()
        params?.delete()
    }

    def 'params should override params-file and override params in the config file' () {
        setup:
        def params = Files.createTempFile('test', '.yml')
        params.text = '''
            alpha: "Hello" 
            beta: "World" 
            omega: "Last"
            '''.stripIndent()
        and:
        def file = Files.createTempFile('test',null)
        file.text = '''
        params {
          alpha = 'x'
        }
        params.beta = 'y'
        params.delta = 'Foo'
        params.gamma = "I'm gamma"
        params.omega = "I'm the last"
        
        process {
          publishDir = [path: params.alpha]
        }
        '''
        when:
        def opt = new CliOptions()
        def run = new CmdRun(paramsFile: params, params: [alpha: 'Hola', beta: 'Mundo'])
        def result = new ConfigBuilder().setOptions(opt).setCmdRun(run).buildGivenFiles(file)

        then:
        result.params.alpha == 'Hola'   // <-- this comes from the CLI
        result.params.beta == 'Mundo'   // <-- this comes from the CLI as well
        result.params.omega == 'Last'   // <-- this comes from the params-file
        result.params.gamma == "I'm gamma"   // <-- from the config
        result.params.delta == 'Foo'         // <-- from the config
        result.process.publishDir == [path: 'Hola']

        cleanup:
        file?.delete()
        params?.delete()
    }

    def 'valid config files' () {

        given:
        def folder = Files.createTempDirectory('test')
        def f1 = folder.resolve('file1')
        def f2 = folder.resolve('file2')

        when:
        new ConfigBuilder()
                    .validateConfigFiles([f1.toString(), f2.toString()])
        then:
        thrown(AbortOperationException)


        when:
        f1.text = '1'; f2.text = '2'
        def files = new ConfigBuilder()
                    .validateConfigFiles([f1.toString(), f2.toString()])
        then:
        files == [f1, f2]


        when:
        files = new ConfigBuilder(homeDir: folder, currentDir: folder)
                    .setUserConfigFiles(f1,f2)
                    .validateConfigFiles()
        then:
        files == [f1, f2]

        cleanup:
        folder.deleteDir()

    }

    def 'should discover default config files' () {
        given:
        def homeDir = Files.createTempDirectory('home')
        def baseDir = Files.createTempDirectory('work')
        def workDir = Files.createTempDirectory('work')

        when:
        def homeConfig = homeDir.resolve('config')
        homeConfig.text = 'foo=1'
        def files1 = new ConfigBuilder(homeDir: homeDir, baseDir: workDir, currentDir: workDir).validateConfigFiles()
        then:
        files1 == [homeConfig]

        when:
        def workConfig = workDir.resolve('nextflow.config')
        workConfig.text = 'bar=2'
        def files2 = new ConfigBuilder(homeDir: homeDir, baseDir: workDir, currentDir: workDir).validateConfigFiles()
        then:
        files2 == [homeConfig, workConfig]

        when:
        def baseConfig = baseDir.resolve('nextflow.config')
        baseConfig.text = 'ciao=3'
        def files3 = new ConfigBuilder(homeDir: homeDir, baseDir: baseDir, currentDir: workDir).validateConfigFiles()
        then:
        files3 == [homeConfig, baseConfig, workConfig]


        cleanup:
        homeDir?.deleteDir()
        workDir?.deleteDir()
    }

    def 'command executor options'() {

        when:
        def opt = new CliOptions()
        def run = new CmdRun(executorOptions: [ alpha: 1, 'beta.x': 'hola', 'beta.y': 'ciao' ])
        def result = new ConfigBuilder().setOptions(opt).setCmdRun(run).buildGivenFiles()
        then:
        result.executor.alpha == 1
        result.executor.beta.x == 'hola'
        result.executor.beta.y == 'ciao'

    }

    def 'run command cluster options'() {

        when:
        def opt = new CliOptions()
        def run = new CmdRun(clusterOptions: [ alpha: 1, 'beta.x': 'hola', 'beta.y': 'ciao' ])
        def result = new ConfigBuilder().setOptions(opt).setCmdRun(run).buildGivenFiles()
        then:
        result.cluster.alpha == 1
        result.cluster.beta.x == 'hola'
        result.cluster.beta.y == 'ciao'

    }

    def 'run with docker'() {

        when:
        def opt = new CliOptions()
        def run = new CmdRun(withDocker: 'cbcrg/piper')
        def config = new ConfigBuilder().setOptions(opt).setCmdRun(run).build()

        then:
        config.docker.enabled
        config.process.container == 'cbcrg/piper'

    }

    def 'run with docker 2'() {

        given:
        def file = Files.createTempFile('test','config')
        file.deleteOnExit()
        file.text =
                '''
                docker {
                    image = 'busybox'
                    enabled = false
                }
                '''

        when:
        def opt = new CliOptions(config: [file.toFile().canonicalPath] )
        def run = new CmdRun(withDocker: '-')
        def config = new ConfigBuilder().setOptions(opt).setCmdRun(run).build()
        then:
        config.docker.enabled
        config.docker.image == 'busybox'
        config.process.container == 'busybox'

        when:
        opt = new CliOptions(config: [file.toFile().canonicalPath] )
        run = new CmdRun(withDocker: 'cbcrg/mybox')
        config = new ConfigBuilder().setOptions(opt).setCmdRun(run).build()
        then:
        config.docker.enabled
        config.process.container == 'cbcrg/mybox'

    }

    def 'run with docker 3'() {
        given:
        def file = Files.createTempFile('test','config')
        file.deleteOnExit()

        when:
        file.text =
                '''
                process.'withName:test'.container = 'busybox'
                '''
        def opt = new CliOptions(config: [file.toFile().canonicalPath])
        def run = new CmdRun(withDocker: '-')
        def config = new ConfigBuilder().setOptions(opt).setCmdRun(run).build()
        then:
        config.docker.enabled
        config.process.'withName:test'.container == 'busybox'

        when:
        file.text =
                '''
                process.container = 'busybox'
                '''
        opt = new CliOptions(config: [file.toFile().canonicalPath])
        run = new CmdRun(withDocker: '-')
        config = new ConfigBuilder().setOptions(opt).setCmdRun(run).build()
        then:
        config.docker.enabled
        config.process.container == 'busybox'

        when:
        opt = new CliOptions()
        run = new CmdRun(withDocker: '-')
        new ConfigBuilder().setOptions(opt).setCmdRun(run).build()
        then:
        def e = thrown(AbortOperationException)
        e.message == 'You have requested to run with Docker but no image was specified'

        when:
        file.text =
                '''
                process.'withName:test'.tag = 'tag'
                '''
        opt = new CliOptions(config: [file.toFile().canonicalPath])
        run = new CmdRun(withDocker: '-')
        new ConfigBuilder().setOptions(opt).setCmdRun(run).build()
        then:
        e = thrown(AbortOperationException)
        e.message == 'You have requested to run with Docker but no image was specified'

    }

    def 'run without docker'() {

        given:
        def file = Files.createTempFile('test','config')
        file.deleteOnExit()
        file.text =
                '''
                docker {
                    image = 'busybox'
                    enabled = true
                }
                '''

        when:
        def opt = new CliOptions(config: [file.toFile().canonicalPath] )
        def run = new CmdRun(withoutDocker: true)
        def config = new ConfigBuilder().setOptions(opt).setCmdRun(run).build()
        then:
        !config.docker.enabled
        config.docker.image == 'busybox'
        !config.process.container

    }

    def 'config with cluster options'() {

        when:
        def opt = new CliOptions()
        def cmd = new CmdNode(clusterOptions: [join: 'x', group: 'y', interface: 'abc', slots: 10, 'tcp.alpha':'uno', 'tcp.beta': 'due'])

        def config = new ConfigBuilder()
                .setOptions(opt)
                .setCmdNode(cmd)
                .build()

        then:
        config.cluster.join == 'x'
        config.cluster.group == 'y'
        config.cluster.interface == 'abc'
        config.cluster.slots == 10
        config.cluster.tcp.alpha == 'uno'
        config.cluster.tcp.beta == 'due'

    }

    def 'has container directive' () {
        when:
        def config = new ConfigBuilder()

        then:
        !config.hasContainerDirective(null)
        !config.hasContainerDirective([:])
        !config.hasContainerDirective([foo: true])
        config.hasContainerDirective([container: 'hello/world'])
        !config.hasContainerDirective([foo: 1, bar: 2])
        !config.hasContainerDirective([foo: 1, bar: 2, baz: [container: 'user/repo']])
        config.hasContainerDirective([foo: 1, bar: 2, $baz: [container: 'user/repo']])
        config.hasContainerDirective([foo: 1, bar: 2, 'withName:baz': [container: 'user/repo']])

    }

    def 'should set session trace options' () {

        given:
        def env = [:]
        def builder = [:] as ConfigBuilder

        when:
        def config = new ConfigObject()
        builder.configRunOptions(config, env, new CmdRun())

        then:
        config.trace instanceof Map
        !config.trace.enabled
        !config.trace.file

        when:
        config = new ConfigObject()
        builder.configRunOptions(config, env, new CmdRun(withTrace: 'some-file'))
        then:
        config.trace instanceof Map
        config.trace.enabled
        config.trace.file == 'some-file'

        when:
        config = new ConfigObject()
        config.trace.file = 'foo.txt'
        builder.configRunOptions(config, env, new CmdRun(withTrace: 'bar.txt'))
        then: // command line should override the config file
        config.trace instanceof Map
        config.trace.enabled
        config.trace.file == 'bar.txt'

        when:
        config = new ConfigObject()
        config.trace.file = 'foo.txt'
        builder.configRunOptions(config, env, new CmdRun(withTrace: '-'))
        then: // command line should override the config file
        config.trace instanceof Map
        config.trace.enabled
        config.trace.file == 'foo.txt'

        when:
        config = new ConfigObject()
        builder.configRunOptions(config, env, new CmdRun(withTrace: '-'))
        then: // command line should override the config file
        config.trace instanceof Map
        config.trace.enabled
        config.trace.file == 'trace-20221001.txt'
    }

    def 'should set session report options' () {

        given:
        def env = [:]
        def builder = [:] as ConfigBuilder

        when:
        def config = new ConfigObject()
        builder.configRunOptions(config, env, new CmdRun())
        then:
        !config.report

        when:
        config = new ConfigObject()
        config.report.file = 'foo.html'
        builder.configRunOptions(config, env, new CmdRun())
        then:
        config.report instanceof Map
        !config.report.enabled
        config.report.file == 'foo.html'

        when:
        config = new ConfigObject()
        builder.configRunOptions(config, env, new CmdRun(withReport: 'my-report.html'))
        then:
        config.report instanceof Map
        config.report.enabled
        config.report.file == 'my-report.html'

        when:
        config = new ConfigObject()
        config.report.file = 'this-report.html'
        builder.configRunOptions(config, env, new CmdRun(withReport: 'my-report.html'))
        then:
        config.report instanceof Map
        config.report.enabled
        config.report.file == 'my-report.html'

        when:
        config = new ConfigObject()
        config.report.file = 'this-report.html'
        builder.configRunOptions(config, env, new CmdRun(withReport: '-'))
        then:
        config.report instanceof Map
        config.report.enabled
        config.report.file == 'this-report.html'

        when:
        config = new ConfigObject()
        builder.configRunOptions(config, env, new CmdRun(withReport: '-'))
        then:
        config.report instanceof Map
        config.report.enabled
        config.report.file == 'report-20221001.html'
    }


    def 'should set session dag options' () {

        given:
        def env = [:]
        def builder = [:] as ConfigBuilder

        when:
        def config = new ConfigObject()
        builder.configRunOptions(config, env, new CmdRun())
        then:
        !config.dag

        when:
        config = new ConfigObject()
        config.dag.file = 'foo-dag.html'
        builder.configRunOptions(config, env, new CmdRun())
        then:
        config.dag instanceof Map
        !config.dag.enabled
        config.dag.file == 'foo-dag.html'

        when:
        config = new ConfigObject()
        builder.configRunOptions(config, env, new CmdRun(withDag: 'my-dag.html'))
        then:
        config.dag instanceof Map
        config.dag.enabled
        config.dag.file == 'my-dag.html'

        when:
        config = new ConfigObject()
        config.dag.file = 'this-dag.html'
        builder.configRunOptions(config, env, new CmdRun(withDag: 'my-dag.html'))
        then:
        config.dag instanceof Map
        config.dag.enabled
        config.dag.file == 'my-dag.html'

        when:
        config = new ConfigObject()
        config.dag.file = 'this-dag.html'
        builder.configRunOptions(config, env, new CmdRun(withDag: '-'))
        then:
        config.dag instanceof Map
        config.dag.enabled
        config.dag.file == 'this-dag.html'

        when:
        config = new ConfigObject()
        builder.configRunOptions(config, env, new CmdRun(withDag: '-'))
        then:
        config.dag instanceof Map
        config.dag.enabled
        config.dag.file == 'dag-20221001.html'
    }

    def 'should set session weblog options' () {

        given:
        def env = [:]
        def builder = [:] as ConfigBuilder

        when:
        def config = new ConfigObject()
        builder.configRunOptions(config, env, new CmdRun())
        then:
        !config.weblog

        when:
        config = new ConfigObject()
        config.weblog.url = 'http://bar.com'
        builder.configRunOptions(config, env, new CmdRun())
        then:
        config.weblog instanceof Map
        !config.weblog.enabled
        config.weblog.url == 'http://bar.com'

        when:
        config = new ConfigObject()
        builder.configRunOptions(config, env, new CmdRun(withWebLog: 'http://foo.com'))
        then:
        config.weblog instanceof Map
        config.weblog.enabled
        config.weblog.url == 'http://foo.com'

        when:
        config = new ConfigObject()
        config.weblog.enabled = true
        config.weblog.url = 'http://bar.com'
        builder.configRunOptions(config, env, new CmdRun(withWebLog: 'http://foo.com'))
        then:
        config.weblog instanceof Map
        config.weblog.enabled
        config.weblog.url == 'http://foo.com'

        when:
        config = new ConfigObject()
        config.weblog.enabled = true
        config.weblog.url = 'http://bar.com'
        builder.configRunOptions(config, env, new CmdRun(withWebLog: '-'))
        then:
        config.weblog instanceof Map
        config.weblog.enabled
        config.weblog.url == 'http://bar.com'

        when:
        config = new ConfigObject()
        config.weblog.enabled = true
        builder.configRunOptions(config, env, new CmdRun(withWebLog: '-'))
        then:
        config.weblog instanceof Map
        config.weblog.enabled
        config.weblog.url == 'http://localhost'

    }

    def 'should set session timeline options' () {

        given:
        def env = [:]
        def builder = [:] as ConfigBuilder

        when:
        def config = new ConfigObject()
        builder.configRunOptions(config, env, new CmdRun())
        then:
        !config.timeline

        when:
        config = new ConfigObject()
        config.timeline.file = 'my-file.html'
        builder.configRunOptions(config, env, new CmdRun())
        then:
        config.timeline instanceof Map
        !config.timeline.enabled
        config.timeline.file == 'my-file.html'

        when:
        config = new ConfigObject()
        builder.configRunOptions(config, env, new CmdRun(withTimeline: 'my-timeline.html'))
        then:
        config.timeline instanceof Map
        config.timeline.enabled
        config.timeline.file == 'my-timeline.html'

        when:
        config = new ConfigObject()
        config.timeline.enabled = true
        config.timeline.file = 'this-timeline.html'
        builder.configRunOptions(config, env, new CmdRun(withTimeline: 'my-timeline.html'))
        then:
        config.timeline instanceof Map
        config.timeline.enabled
        config.timeline.file == 'my-timeline.html'

        when:
        config = new ConfigObject()
        config.timeline.enabled = true
        config.timeline.file = 'my-timeline.html'
        builder.configRunOptions(config, env, new CmdRun(withTimeline: '-'))
        then:
        config.timeline instanceof Map
        config.timeline.enabled
        config.timeline.file == 'my-timeline.html'

        when:
        config = new ConfigObject()
        config.timeline.enabled = true
        builder.configRunOptions(config, env, new CmdRun(withTimeline: '-'))
        then:
        config.timeline instanceof Map
        config.timeline.enabled
        config.timeline.file == 'timeline-20221001.html'
    }

    def 'should set tower options' () {

        given:
        def env = [:]
        def builder = [:] as ConfigBuilder

        when:
        def config = new ConfigObject()
        builder.configRunOptions(config, env, new CmdRun())
        then:
        !config.tower

        when:
        config = new ConfigObject()
        config.tower.endpoint = 'http://foo.com'
        builder.configRunOptions(config, env, new CmdRun())
        then:
        config.tower instanceof Map
        !config.tower.enabled
        config.tower.endpoint == 'http://foo.com'

        when:
        config = new ConfigObject()
        builder.configRunOptions(config, env, new CmdRun(withTower: 'http://bar.com'))
        then:
        config.tower instanceof Map
        config.tower.enabled
        config.tower.endpoint == 'http://bar.com'

        when:
        config = new ConfigObject()
        config.tower.endpoint = 'http://foo.com'
        builder.configRunOptions(config, env, new CmdRun(withTower: '-'))
        then:
        config.tower instanceof Map
        config.tower.enabled
        config.tower.endpoint == 'http://foo.com'

        when:
        config = new ConfigObject()
        builder.configRunOptions(config, env, new CmdRun(withTower: '-'))
        then:
        config.tower instanceof Map
        config.tower.enabled
        config.tower.endpoint == 'https://api.cloud.seqera.io'
    }

    def 'should set wave options' () {

        given:
        def env = [:]
        def builder = [:] as ConfigBuilder

        when:
        def config = new ConfigObject()
        builder.configRunOptions(config, env, new CmdRun())
        then:
        !config.wave

        when:
        config = new ConfigObject()
        config.wave.endpoint = 'http://foo.com'
        builder.configRunOptions(config, env, new CmdRun())
        then:
        config.wave instanceof Map
        !config.wave.enabled
        config.wave.endpoint == 'http://foo.com'

        when:
        config = new ConfigObject()
        builder.configRunOptions(config, env, new CmdRun(withWave: 'http://bar.com'))
        then:
        config.wave instanceof Map
        config.wave.enabled
        config.wave.endpoint == 'http://bar.com'

        when:
        config = new ConfigObject()
        config.wave.endpoint = 'http://foo.com'
        builder.configRunOptions(config, env, new CmdRun(withWave: '-'))
        then:
        config.wave instanceof Map
        config.wave.enabled
        config.wave.endpoint == 'http://foo.com'

        when:
        config = new ConfigObject()
        builder.configRunOptions(config, env, new CmdRun(withWave: '-'))
        then:
        config.wave instanceof Map
        config.wave.enabled
        config.wave.endpoint == 'https://wave.seqera.io'
    }

    def 'should set cloudcache options' () {

        given:
        def env = [:]
        def builder = [:] as ConfigBuilder

        when:
        def config = new ConfigObject()
        builder.configRunOptions(config, env, new CmdRun())
        then:
        !config.cloudcache

        when:
        config = new ConfigObject()
        config.cloudcache.path = 's3://foo/bar'
        builder.configRunOptions(config, env, new CmdRun())
        then:
        config.cloudcache instanceof Map
        !config.cloudcache.enabled
        config.cloudcache.path == 's3://foo/bar'

        when:
        config = new ConfigObject()
        builder.configRunOptions(config, env, new CmdRun(cloudCachePath: 's3://this/that'))
        then:
        config.cloudcache instanceof Map
        config.cloudcache.enabled
        config.cloudcache.path == 's3://this/that'

        when:
        config = new ConfigObject()
        builder.configRunOptions(config, env, new CmdRun(cloudCachePath: '-'))
        then:
        config.cloudcache instanceof Map
        config.cloudcache.enabled
        !config.cloudcache.path

        when:
        config = new ConfigObject()
        config.cloudcache.path = 's3://alpha/delta'
        builder.configRunOptions(config, env, new CmdRun(cloudCachePath: '-'))
        then:
        config.cloudcache instanceof Map
        config.cloudcache.enabled
        config.cloudcache.path == 's3://alpha/delta'

        when:
        config = new ConfigObject()
        config.cloudcache.path = 's3://alpha/delta'
        builder.configRunOptions(config, env, new CmdRun(cloudCachePath: 's3://should/override/config'))
        then:
        config.cloudcache instanceof Map
        config.cloudcache.enabled
        config.cloudcache.path == 's3://should/override/config'

        when:
        config = new ConfigObject()
        config.cloudcache.enabled = false
        builder.configRunOptions(config, env, new CmdRun(cloudCachePath: 's3://should/override/config'))
        then:
        config.cloudcache instanceof Map
        !config.cloudcache.enabled
        config.cloudcache.path == 's3://should/override/config'

        when:
        config = new ConfigObject()
        builder.configRunOptions(config, [NXF_CLOUDCACHE_PATH:'s3://foo'], new CmdRun(cloudCachePath: 's3://should/override/env'))
        then:
        config.cloudcache instanceof Map
        config.cloudcache.enabled
        config.cloudcache.path == 's3://should/override/env'

        when:
        config = new ConfigObject()
        config.cloudcache.path = 's3://config/path'
        builder.configRunOptions(config, [NXF_CLOUDCACHE_PATH:'s3://foo'], new CmdRun())
        then:
        config.cloudcache instanceof Map
        config.cloudcache.enabled
        config.cloudcache.path == 's3://config/path'

        when:
        config = new ConfigObject()
        config.cloudcache.path = 's3://config/path'
        builder.configRunOptions(config, [NXF_CLOUDCACHE_PATH:'s3://foo'], new CmdRun(cloudCachePath: 's3://should/override/config'))
        then:
        config.cloudcache instanceof Map
        config.cloudcache.enabled
        config.cloudcache.path == 's3://should/override/config'

    }

    def 'should enable conda env' () {

        given:
        def env = [:]
        def builder = [:] as ConfigBuilder

        when:
        def config = new ConfigObject()
        builder.configRunOptions(config, env, new CmdRun())
        then:
        !config.conda

        when:
        config = new ConfigObject()
        config.conda.createOptions = 'something'
        builder.configRunOptions(config, env, new CmdRun())
        then:
        config.conda instanceof Map
        !config.conda.enabled
        config.conda.createOptions == 'something'

        when:
        config = new ConfigObject()
        builder.configRunOptions(config, env, new CmdRun(withConda: 'my-recipe.yml'))
        then:
        config.conda instanceof Map
        config.conda.enabled
        config.process.conda == 'my-recipe.yml'

        when:
        config = new ConfigObject()
        config.conda.enabled = true
        builder.configRunOptions(config, env, new CmdRun(withConda: 'my-recipe.yml'))
        then:
        config.conda instanceof Map
        config.conda.enabled
        config.process.conda == 'my-recipe.yml'

        when:
        config = new ConfigObject()
        config.process.conda = 'my-recipe.yml'
        builder.configRunOptions(config, env, new CmdRun(withConda: '-'))
        then:
        config.conda instanceof Map
        config.conda.enabled
        config.process.conda == 'my-recipe.yml'
    }

    def 'should disable conda env' () {
        given:
        def file = Files.createTempFile('test','config')
        file.deleteOnExit()
        file.text =
                '''
                conda {
                    enabled = true
                }
                '''

        when:
        def opt = new CliOptions(config: [file.toFile().canonicalPath] )
        def run = new CmdRun(withoutConda: true)
        def config = new ConfigBuilder().setOptions(opt).setCmdRun(run).build()
        then:
        !config.conda.enabled
        !config.process.conda
    }

    def 'should enable spack env' () {

        given:
        def env = [:]
        def builder = [:] as ConfigBuilder

        when:
        def config = new ConfigObject()
        builder.configRunOptions(config, env, new CmdRun())
        then:
        !config.spack

        when:
        config = new ConfigObject()
        config.spack.createOptions = 'something'
        builder.configRunOptions(config, env, new CmdRun())
        then:
        config.spack instanceof Map
        !config.spack.enabled
        config.spack.createOptions == 'something'

        when:
        config = new ConfigObject()
        builder.configRunOptions(config, env, new CmdRun(withSpack: 'my-recipe.yaml'))
        then:
        config.spack instanceof Map
        config.spack.enabled
        config.process.spack == 'my-recipe.yaml'

        when:
        config = new ConfigObject()
        config.spack.enabled = true
        builder.configRunOptions(config, env, new CmdRun(withSpack: 'my-recipe.yaml'))
        then:
        config.spack instanceof Map
        config.spack.enabled
        config.process.spack == 'my-recipe.yaml'

        when:
        config = new ConfigObject()
        config.process.spack = 'my-recipe.yaml'
        builder.configRunOptions(config, env, new CmdRun(withSpack: '-'))
        then:
        config.spack instanceof Map
        config.spack.enabled
        config.process.spack == 'my-recipe.yaml'
    }

    def 'should disable spack env' () {
        given:
        def file = Files.createTempFile('test','config')
        file.deleteOnExit()
        file.text =
                '''
                spack {
                    enabled = true
                }
                '''

        when:
        def opt = new CliOptions(config: [file.toFile().canonicalPath] )
        def run = new CmdRun(withoutSpack: true)
        def config = new ConfigBuilder().setOptions(opt).setCmdRun(run).build()
        then:
        !config.spack.enabled
        !config.process.spack
    }

    def 'SHOULD SET `RESUME` OPTION'() {

        given:
        def env = [:]
        def builder = [:] as ConfigBuilder

        when:
        def config = new ConfigObject()
        builder.configRunOptions(config, env, new CmdRun())
        then:
        !config.isSet('resume')

        when:
        config = new ConfigObject()
        config.resume ='alpha-beta-delta'
        builder.configRunOptions(config, env, new CmdRun())
        then:
        config.resume == 'alpha-beta-delta'

        when:
        config = new ConfigObject()
        config.resume ='alpha-beta-delta'
        builder.configRunOptions(config, env, new CmdRun(resume: 'xxx-yyy'))
        then:
        config.resume == 'xxx-yyy'

        when:
        config = new ConfigObject()
        config.resume ='this-that'
        builder.configRunOptions(config, env, new CmdRun(resume: 'xxx-yyy'))
        then:
        config.resume == 'xxx-yyy'
    }

    def 'should set `workDir`' () {

        given:
        def config = new ConfigObject()
        def builder = [:] as ConfigBuilder

        when:
        builder.configRunOptions(config, [:], new CmdRun())
        then:
        config.workDir == 'work'

        when:
        config = new ConfigObject()
        builder.configRunOptions(config, [NXF_WORK: '/foo/bar'], new CmdRun())
        then:
        config.workDir == '/foo/bar'

        when:
        config = new ConfigObject()
        config.workDir = 'hello/there'
        builder.configRunOptions(config, [:], new CmdRun())
        then:
        config.workDir == 'hello/there'

        when:
        config = new ConfigObject()
        config.workDir = 'hello/there'
        builder.configRunOptions(config, [:], new CmdRun(workDir: 'my/work/dir'))
        then:
        config.workDir == 'my/work/dir'
    }

    def 'should set `libDir`' () {
        given:
        def config = new ConfigObject()
        def builder = [:] as ConfigBuilder

        when:
        builder.configRunOptions(config, [:], new CmdRun())
        then:
        !config.isSet('libDir')

        when:
        builder.configRunOptions(config, [NXF_LIB:'/foo/bar'], new CmdRun())
        then:
        config.libDir == '/foo/bar'

        when:
        builder.configRunOptions(config, [:], new CmdRun(libPath: 'my/lib/dir'))
        then:
        config.libDir == 'my/lib/dir'
    }

    def 'should set `cacheable`' () {
        given:
        def env = [:]
        def config
        def builder = [:] as ConfigBuilder

        when:
        config = new ConfigObject()
        builder.configRunOptions(config, env, new CmdRun())
        then:
        !config.isSet('cacheable')

        when:
        config = new ConfigObject()
        builder.configRunOptions(config, env, new CmdRun(cacheable: false))
        then:
        config.cacheable == false

        when:
        config = new ConfigObject()
        builder.configRunOptions(config, env, new CmdRun(cacheable: true))
        then:
        config.cacheable == true
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

    def 'should set profile options' () {

        def builder

        when:
        builder = new ConfigBuilder().setCmdRun(new CmdRun(profile: 'foo'))
        then:
        builder.profile == 'foo'
        builder.validateProfile

        when:
        builder = new ConfigBuilder().setCmdRun(new CmdRun())
        then:
        builder.profile == 'standard'
        !builder.validateProfile

        when:
        builder = new ConfigBuilder().setCmdRun(new CmdRun(profile: 'standard'))
        then:
        builder.profile == 'standard'
        builder.validateProfile
    }

    def 'should set config options' () {
        def builder

        when:
        builder = new ConfigBuilder().setCmdConfig(new CmdConfig())
        then:
        !builder.showAllProfiles

        when:
        builder = new ConfigBuilder().setCmdConfig(new CmdConfig(showAllProfiles: true))
        then:
        builder.showAllProfiles

        when:
        builder = new ConfigBuilder().setCmdConfig(new CmdConfig(profile: 'foo'))
        then:
        builder.profile == 'foo'
        builder.validateProfile

    }

    def 'should set params into config object' () {

        given:
        def emptyFile = Files.createTempFile('empty','config').toFile()
        def EMPTY = [emptyFile.toString()]

        def configFile = Files.createTempFile('test','config').toFile()
        configFile.deleteOnExit()
        configFile.text = '''
          params.foo = 1
          params.bar = 2
          params.data = '/some/path'
        '''
        configFile = configFile.toString()

        def jsonFile = Files.createTempFile('test','.json').toFile()
        jsonFile.text = '''
          {
            "foo": 10,
            "bar": 20
          }
        '''
        jsonFile = jsonFile.toString()

        def yamlFile = Files.createTempFile('test','.yaml').toFile()
        yamlFile.text = '''
          {
            "foo": 100,
            "bar": 200
          }
        '''
        yamlFile = yamlFile.toString()

        def config

        when:
        config = new ConfigBuilder().setOptions(new CliOptions(config: EMPTY)).setCmdRun(new CmdRun()).build()
        then:
        config.params == [:]

        // get params for the CLI
        when:
        config = new ConfigBuilder().setOptions(new CliOptions(config: EMPTY)).setCmdRun(new CmdRun(params: [foo:'one', bar:'two'])).build()
        then:
        config.params == [foo:'one', bar:'two']

        // get params from config file
        when:
        config = new ConfigBuilder().setOptions(new CliOptions(config: [configFile])).setCmdRun(new CmdRun()).build()
        then:
        config.params == [foo:1, bar:2, data: '/some/path']

        // get params form JSON file
        when:
        config = new ConfigBuilder().setOptions(new CliOptions(config: EMPTY)).setCmdRun(new CmdRun(paramsFile: jsonFile)).build()
        then:
        config.params == [foo:10, bar:20]

        // get params from YAML file
        when:
        config = new ConfigBuilder().setOptions(new CliOptions(config: EMPTY)).setCmdRun(new CmdRun(paramsFile: yamlFile)).build()
        then:
        config.params == [foo:100, bar:200]

        // cli override config
        when:
        config = new ConfigBuilder().setOptions(new CliOptions(config: [configFile])).setCmdRun(new CmdRun(params:[foo:'hello', baz:'world'])).build()
        then:
        config.params == [foo:'hello', bar:2, baz: 'world', data: '/some/path']

        // CLI override JSON
        when:
        config = new ConfigBuilder().setOptions(new CliOptions(config: EMPTY)).setCmdRun(new CmdRun(params:[foo:'hello', baz:'world'], paramsFile: jsonFile)).build()
        then:
        config.params == [foo:'hello', bar:20, baz: 'world']

        // JSON override config
        when:
        config = new ConfigBuilder().setOptions(new CliOptions(config: [configFile])).setCmdRun(new CmdRun(paramsFile: jsonFile)).build()
        then:
        config.params == [foo:10, bar:20, data: '/some/path']


        // CLI override JSON that override config
        when:
        config = new ConfigBuilder().setOptions(new CliOptions(config: [configFile])).setCmdRun(new CmdRun(paramsFile: jsonFile, params: [foo:'Ciao'])).build()
        then:
        config.params == [foo:'Ciao', bar:20, data: '/some/path']
    }

    def 'should run with conda' () {

        when:
        def config = new ConfigBuilder().setCmdRun(new CmdRun(withConda: '/some/path/env.yml')).build()
        then:
        config.process.conda == '/some/path/env.yml'

    }

    def 'should run with spack' () {

        when:
        def config = new ConfigBuilder().setCmdRun(new CmdRun(withSpack: '/some/path/env.yaml')).build()
        then:
        config.process.spack == '/some/path/env.yaml'

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
        def opt = new CliOptions(config: [file.toFile().canonicalPath] )
        def cfg = new ConfigBuilder().setOptions(opt).build()
        then:
        cfg.params.foo == System.getenv('HOME')

        when:
        file.text =
                '''
                params.foo = bar
                '''
        opt = new CliOptions(config: [file.toFile().canonicalPath] )
        new ConfigBuilder().setOptions(opt).build()
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
        def opt = new CliOptions(config: [file.toFile().canonicalPath] )
        def builder = new ConfigBuilder()
                .setOptions(opt)
                .showMissingVariables(true)
        def cfg = builder.buildConfigObject()
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
        def opt = new CliOptions(config: [file.toFile().canonicalPath] )
        new ConfigBuilder().setOptions(opt).build()
        then:
        def e = thrown(ConfigParseException)
        e.message == "Unknown config attribute `foo.bar` -- check config file: ${file.toRealPath()}".toString()

    }


    def 'should collect config files' () {

        given:
        def slurper = new ConfigParser()
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


    def 'should configure notification' () {

        given:
        Map config

        when:
        config = new ConfigBuilder().setCmdRun(new CmdRun()).build()
        then:
        !config.notification

        when:
        config = new ConfigBuilder().setCmdRun(new CmdRun(withNotification: true)).build()
        then:
        config.notification.enabled == true

        when:
        config = new ConfigBuilder().setCmdRun(new CmdRun(withNotification: false)).build()
        then:
        config.notification.enabled == false
        config.notification.to == null

        when:
        config = new ConfigBuilder().setCmdRun(new CmdRun(withNotification: 'yo@nextflow.com')).build()
        then:
        config.notification.enabled == true
        config.notification.to == 'yo@nextflow.com'
    }

    def 'should configure fusion' () {

        given:
        Map config

        when:
        config = new ConfigBuilder().setCmdRun(new CmdRun()).build()
        then:
        !config.fusion

        when:
        config = new ConfigBuilder().setCmdRun(new CmdRun(withFusion: true)).build()
        then:
        config.fusion.enabled == true

        when:
        config = new ConfigBuilder().setCmdRun(new CmdRun(withFusion: false)).build()
        then:
        config.fusion == [enabled: false]

        when:
        config = new ConfigBuilder().setCmdRun(new CmdRun(withFusion: true)).build()
        then:
        config.fusion == [enabled: true]
    }

    def 'should configure stub run mode' () {
        given:
        Map config

        when:
        config = new ConfigBuilder().setCmdRun(new CmdRun()).build()
        then:
        !config.stubRun

        when:
        config = new ConfigBuilder().setCmdRun(new CmdRun(stubRun: true)).build()
        then:
        config.stubRun == true
    }

    def 'should configure preview mode' () {
        given:
        Map config

        when:
        config = new ConfigBuilder().setCmdRun(new CmdRun()).build()
        then:
        !config.preview

        when:
        config = new ConfigBuilder().setCmdRun(new CmdRun(preview: true)).build()
        then:
        config.preview == true
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
        result = builder.buildConfig0(ENV, [CONFIG])
        then:
        result.process.container == 'base'
        result.process.executor == 'sge'
        result.process.queue == 'short'

        when:
        builder.profile = 'cfg2'
        result = builder.buildConfig0(ENV, [CONFIG])
        then:
        result.process.container == 'base'
        result.process.executor == 'batch'
        result.process.queue == 'long'

        when:
        builder.profile = 'cfg1,docker'
        result = builder.buildConfig0(ENV, [CONFIG])
        then:
        result.process.container ==  'foo/1'
        result.process.executor == 'sge'
        result.process.queue == 'short'
        result.docker.enabled
        !result.isSet('singularity')

        when:
        builder.profile = 'cfg1,singularity'
        result = builder.buildConfig0(ENV, [CONFIG])
        then:
        result.process.container ==  'bar-2.img'
        result.process.executor == 'sge'
        result.process.queue == 'cn-el7'
        result.singularity.enabled
        !result.isSet('docker')

        when:
        builder.profile = 'cfg2,singularity'
        result = builder.buildConfig0(ENV, [CONFIG])
        then:
        result.process.container ==  'bar-2.img'
        result.process.executor == 'batch'
        result.process.queue == 'cn-el7'
        result.singularity.enabled
        !result.isSet('docker')

        when:
        builder.profile = 'missing'
        builder.buildConfig0(ENV, [CONFIG])
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
        result = builder.buildConfig0(ENV, [CONFIG])
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
        def cfg1 = new ConfigBuilder().buildConfig0([:], [file1])
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

        def cfg2 = new ConfigBuilder().buildConfig0([:], [file2])
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
        def cfg1 = new ConfigBuilder().buildConfig0([:], [file1])
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
        def cfg1 = new ConfigBuilder().buildConfig0([:], [file1])
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
        def cfg = new ConfigBuilder().setProfile('foo').buildConfig0([:], [file1])
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
        def cfg = new ConfigBuilder().setProfile('foo').buildConfig0([:], [file1])
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
        def cfg = new ConfigBuilder().setProfile('foo,bar') .buildConfig0([:], [file1])
        then:
        cfg.process.cpus == 1
        cfg.params.alpha == 1
        cfg.params.delta == 20
        cfg.params.gamma == 30

        cleanup:
        folder?.deleteDir()
    }

    def 'CLI params should overwrite only the key provided when nested'() {
        given:
        def folder = File.createTempDir()
        def configMain = new File(folder, 'nextflow.config').absoluteFile

        configMain.text = """
        params {
            foo = 'Hello'
            bar = "Monde"
            baz {
                x = "Ciao"
                y = "mundo"
                z { 
                    alpha = "Hallo"
                    beta  = "World"
                }
            }
            
        }        
        """

        when:
        def opt = new CliOptions()
        def run = new CmdRun(params: [bar: "world", 'baz.y': "mondo", 'baz.z.beta': "Welt"])
        def config = new ConfigBuilder(env: [NXF_CONFIG_FILE: configMain.toString()]).setOptions(opt).setCmdRun(run).build()

        then:
        config.params.foo == 'Hello'
        config.params.bar == 'world'
        config.params.baz.x == 'Ciao'
        config.params.baz.y == 'mondo'
        //tests recursion
        config.params.baz.z.alpha == 'Hallo'
        config.params.baz.z.beta == 'Welt'

        cleanup:
        folder?.deleteDir()
    }

    @Unroll
    def 'should merge config params' () {
        given:
        def builder = new ConfigBuilder()

        expect:
        def cfg = new ConfigObject(); if(CONFIG) cfg.putAll(CONFIG)
        and:
        builder.mergeMaps(cfg, PARAMS, false) == EXPECTED

        where:
        CONFIG                      | PARAMS | EXPECTED
        [foo:1]                     | null                          | [foo:1]
        null                        | [bar:2]                       | [bar:2]
        [foo:1]                     | [bar:2]                       | [foo: 1, bar: 2]
        [foo:1]                     | [bar:null]                    | [foo: 1, bar: null]
        [foo:1]                     | [foo:null]                    | [foo: null]
        [foo:1, bar:[:]]            | [bar: [x:10, y:20]]           | [foo: 1, bar: [x:10, y:20]]
        [foo:1, bar:[x:1, y:2]]     | [bar: [x:10, y:20]]           | [foo: 1, bar: [x:10, y:20]]
        [foo:1, bar:[x:1, y:2]]     | [foo: 2, bar: [x:10, y:20]]   | [foo: 2, bar: [x:10, y:20]]
        [foo:1, bar:null]           | [bar: [x:10, y:20]]           | [foo: 1, bar: [x:10, y:20]]
        [foo:1, bar:2]              | [bar: [x:10, y:20]]           | [foo: 1, bar: [x:10, y:20]]
        [foo:1, bar:[x:1, y:2]]     | [bar: 2]                      | [foo: 1, bar: 2]
    }

    @Unroll
    def 'should merge config strict params' () {
        given:
        def builder = new ConfigBuilder()

        expect:
        def cfg = new ConfigObject(); if(CONFIG) cfg.putAll(CONFIG)
        and:
        builder.mergeMaps(cfg, PARAMS, true) == EXPECTED

        where:
        CONFIG                      | PARAMS                        | EXPECTED
        [:]                         | [bar:2]                       | [bar:2]
        [foo:1]                     | null                          | [foo:1]
        null                        | [bar:2]                       | [bar:2]
        [foo:1]                     | [bar:2]                       | [foo: 1, bar: 2]
        [foo:1]                     | [bar:null]                    | [foo: 1, bar: null]
        [foo:1]                     | [foo:null]                    | [foo: null]
        [foo:1, bar:[:]]            | [bar: [x:10, y:20]]           | [foo: 1, bar: [x:10, y:20]]
        [foo:1, bar:[x:1, y:2]]     | [bar: [x:10, y:20]]           | [foo: 1, bar: [x:10, y:20]]
        [foo:1, bar:[x:1, y:2]]     | [foo: 2, bar: [x:10, y:20]]   | [foo: 2, bar: [x:10, y:20]]
    }

    def 'prevent config side effects' () {
        given:
        def folder = Files.createTempDirectory('test')
        and:
        def config = folder.resolve('nf.config')
        config.text = '''\
        params.test.foo = "foo_def"
        params.test.bar = "bar_def"        
        '''.stripIndent()

        when:
        def cfg1 = new ConfigBuilder()
                .setOptions( new CliOptions(userConfig: [config.toString()]))
                .build()
        then:
        cfg1.params.test.foo == "foo_def"
        cfg1.params.test.bar == "bar_def"

        
        when:
        def cfg2 = new ConfigBuilder()
                .setOptions( new CliOptions(userConfig: [config.toString()]))
                .setCmdRun( new CmdRun(params: ['test.foo': 'CLI_FOO'] ))
                .build()
        then:
        cfg2.params.test.foo == "CLI_FOO"
        cfg2.params.test.bar == "bar_def"

        cleanup:
        folder?.deleteDir()
    }

    def 'parse nested json' () {
        given:
        def folder = Files.createTempDirectory('test')
        and:
        def config = folder.resolve('nf.json')
        config.text = '''\
        {
            "title": "something",
            "nested": {
                "name": "Mike",
                "and": {
                    "more": "nesting",
                    "still": {
                        "another": "layer"
                    }
                }
            }
        }
        '''.stripIndent()

        when:
        def cfg1 = new ConfigBuilder().setCmdRun(new CmdRun(paramsFile: config.toString())).build()

        then:
        cfg1.params.title == "something"
        cfg1.params.nested.name == 'Mike'
        cfg1.params.nested.and.more == 'nesting'
        cfg1.params.nested.and.still.another == 'layer'

        cleanup:
        folder?.deleteDir()
    }

    def 'parse nested yaml' () {
        given:
        def folder = Files.createTempDirectory('test')
        and:
        def config = folder.resolve('nf.yaml')
        config.text = '''\
            title: "something"
            nested: 
              name: "Mike"
              and:
                more: nesting
                still:
                  another: layer      
        '''.stripIndent()

        when:
        def cfg1 = new ConfigBuilder().setCmdRun(new CmdRun(paramsFile: config.toString())).build()

        then:
        cfg1.params.title == "something"
        cfg1.params.nested.name == 'Mike'
        cfg1.params.nested.and.more == 'nesting'
        cfg1.params.nested.and.still.another == 'layer'

        cleanup:
        folder?.deleteDir()
    }


    def 'should return parsed config' () {
        given:
        def cmd = new CmdRun(profile: 'first', withTower: 'http://foo.com', launcher: new Launcher())
        def base = Files.createTempDirectory('test')
        base.resolve('nextflow.config').text = '''
        profiles {
            first {
                params {
                  foo = 'Hello world'
                  awsKey = 'xyz'
                }
                process {
                    executor = { 'local' }
                }
            }
            second {
                params.none = 'Blah'
            }
        }
        '''
        when:
        def txt = ConfigBuilder.resolveConfig(base, cmd)
        then:
        txt == '''\
            params {
               foo = 'Hello world'
               awsKey = '[secret]'
            }
            
            process {
               executor = { 'local' }
            }

            workDir = 'work'
            
            tower {
               enabled = true
               endpoint = 'http://foo.com'
            }
            '''.stripIndent()

        cleanup:
        base?.deleteDir()
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
        def cfg = new ConfigBuilder().setProfile('test').buildConfig0([:], [main])
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
        def cfg = new ConfigBuilder().setBaseDir(folder).buildConfig0([:], [text])
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
        def opt = new CliOptions()
        def config = new ConfigBuilder().setBaseDir(folder).setOptions(opt).buildGivenFiles(configMain)
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
        def opt = new CliOptions()
        def config = new ConfigBuilder().setBaseDir(folder).setOptions(opt).buildGivenFiles(configMain)
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
        def opt = new CliOptions()
        def config = new ConfigBuilder()
                            .setBaseDir(folder)
                            .setOptions(opt)
                            .setStripSecrets(true)
                            .buildGivenFiles(configMain)
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

