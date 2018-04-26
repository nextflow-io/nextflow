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

package nextflow.config
import java.nio.file.Files
import java.nio.file.Paths

import nextflow.cli.CliOptions
import nextflow.cli.CmdConfig
import nextflow.cli.CmdNode
import nextflow.cli.CmdRun
import nextflow.exception.AbortOperationException
import nextflow.exception.ConfigParseException
import spock.lang.Specification
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
        }
        '''

        when:
        def cfg = builder.buildConfig0([:], [text])
        then:
        cfg.params.p == '/base/path/1'
        cfg.params.q == '/base/path/2'

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
        def result = new ConfigBuilder().setOptions(opt).setCmdRun(run).buildConfig([file])

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
        def config = new ConfigBuilder().setOptions(opt).setCmdRun(run).buildConfig([configMain.toPath()])

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
        def config = new ConfigBuilder().setOptions(opt).setCmdRun(run).buildConfig([file])
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
        config = new ConfigBuilder().setOptions(opt).setCmdRun(run).buildConfig([file])
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
        config = new ConfigBuilder().setOptions(opt).setCmdRun(run).buildConfig([file])
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

    def 'valid config files' () {

        given:
        def path = Files.createTempDirectory('test')

        when:
        def f1 = path.resolve('file1')
        def f2 = path.resolve('file2')
        def files = new ConfigBuilder().validateConfigFiles([f1.toString(), f2.toString()])
        then:
        thrown(AbortOperationException)

        when:
        f1.text = '1'; f2.text = '2'
        files = new ConfigBuilder().validateConfigFiles([f1.toString(), f2.toString()])
        then:
        files == [f1, f2]

        cleanup:
        path.deleteDir()

    }

    def 'command executor options'() {

        when:
        def opt = new CliOptions()
        def run = new CmdRun(executorOptions: [ alpha: 1, 'beta.x': 'hola', 'beta.y': 'ciao' ])
        def result = new ConfigBuilder().setOptions(opt).setCmdRun(run).buildConfig([])
        then:
        result.executor.alpha == 1
        result.executor.beta.x == 'hola'
        result.executor.beta.y == 'ciao'

    }

    def 'run command cluster options'() {

        when:
        def opt = new CliOptions()
        def run = new CmdRun(clusterOptions: [ alpha: 1, 'beta.x': 'hola', 'beta.y': 'ciao' ])
        def result = new ConfigBuilder().setOptions(opt).setCmdRun(run).buildConfig([])
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
                process.$test.container = 'busybox'
                '''
        def opt = new CliOptions(config: [file.toFile().canonicalPath])
        def run = new CmdRun(withDocker: '-')
        def config = new ConfigBuilder().setOptions(opt).setCmdRun(run).build()
        then:
        config.docker.enabled
        config.process.$test.container == 'busybox'

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
        e.message == 'You have requested to run with Docker but no image were specified'

        when:
        file.text =
                '''
                process.$test.tag = 'tag'
                '''
        opt = new CliOptions(config: [file.toFile().canonicalPath])
        run = new CmdRun(withDocker: '-')
        new ConfigBuilder().setOptions(opt).setCmdRun(run).build()
        then:
        e = thrown(AbortOperationException)
        e.message == 'You have requested to run with Docker but no image were specified'

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

    }

    def 'should set session trace options' () {

        given:
        def env = [:]
        def config = new ConfigObject()
        def builder = [:] as ConfigBuilder

        when:
        config.trace.enabled = true
        builder.configRunOptions(config, env, new CmdRun())

        then:
        config.trace instanceof Map
        config.trace.enabled
        !config.trace.file

        when:
        builder.configRunOptions(config, env, new CmdRun(withTrace: 'some-file'))
        then:
        config.trace instanceof Map
        config.trace.enabled
        config.trace.file == 'some-file'

    }

    def 'should set session timeline options' () {

        given:
        def env = [:]
        def config = new ConfigObject()
        def builder = [:] as ConfigBuilder

        when:
        builder.configRunOptions(config, env, new CmdRun())
        then:
        !config.timeline

        when:
        config.timeline.enabled = true
        builder.configRunOptions(config, env, new CmdRun())

        then:
        config.timeline instanceof Map
        config.timeline.enabled
        !config.timeline.file

        when:
        builder.configRunOptions(config, env, new CmdRun(withTimeline: 'my-timeline.html'))
        then:
        config.timeline instanceof Map
        config.timeline.enabled
        config.timeline.file == 'my-timeline.html'

    }


    def 'SHOULD SET `RESUME` OPTION'() {

        given:
        def env = [:]
        def config = new ConfigObject()
        def builder = [:] as ConfigBuilder

        when:
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
        config.libDir == null

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
        config.cacheable == true

        when:
        config = new ConfigObject()
        builder.configRunOptions(config, env, new CmdRun(cacheable: false))
        then:
        config.cacheable == false
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
        e.message == "Unknown config attribute: bar -- check config file: ${file.toRealPath()}".toString()

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
        builder.configFiles == [file1, file2]

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

}

