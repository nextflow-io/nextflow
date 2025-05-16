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

package nextflow.config.parser.v2

import java.nio.file.Files
import java.nio.file.NoSuchFileException
import java.nio.file.Path

import com.sun.net.httpserver.Headers
import com.sun.net.httpserver.HttpExchange
import com.sun.net.httpserver.HttpHandler
import com.sun.net.httpserver.HttpServer
import nextflow.SysEnv
import nextflow.config.ConfigBuilder
import nextflow.config.ConfigClosurePlaceholder
import nextflow.exception.ConfigParseException
import nextflow.util.Duration
import nextflow.util.MemoryUnit
import spock.lang.Specification

import static test.TestHelper.createInMemTempFile

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ConfigParserV2Test extends Specification {

    def 'should get an environment variable' () {
        given:
        SysEnv.push(MAX_CPUS: '1')

        when:
        def CONFIG = '''
            process.cpus = env('MAX_CPUS')
            '''
        def config = new ConfigParserV2().parse(CONFIG)

        then:
        config.process.cpus == '1'

        cleanup:
        SysEnv.pop()
    }

    def 'should parse plugin ids' () {
        given:
        def CONFIG = '''
            plugins {
                id 'foo'
                id 'bar'
                id 'bar'
            }

            process {
                cpus = 1
                mem = 2
            }
            '''

        when:
        def config = new ConfigParserV2().parse(CONFIG)

        then:
        config.plugins == ['foo','bar'] as Set
        and:
        config.process.cpus == 1
    }

    def 'should parse composed config files' () {

        given:
        def folder = Files.createTempDirectory('test')
        def main = folder.resolve('nextflow.config')
        def snippet1 = folder.resolve('snippet1.config')
        def snippet2 = folder.resolve('snippet2.config')

        main.text = """
            profiles {
                alpha {
                    process.disk = '1TB'
                    includeConfig "$snippet1"
                }
            }
            """

        snippet1.text = """
            process.cpus = 4
            process.memory = '8GB'
            includeConfig "$snippet2"
            """

        snippet2.text = '''
            process.ext.foo = 1
            process.ext.bar = 2
            '''

        when:
        def config = new ConfigParserV2().parse(main)
        then:
        config.profiles.alpha.process.cpus == 4
        config.profiles.alpha.process.memory == '8GB'
        config.profiles.alpha.process.disk == '1TB'
        config.profiles.alpha.process.ext.foo == 1
        config.profiles.alpha.process.ext.bar == 2

        cleanup:
        folder?.deleteDir()

    }

    def 'should parse include config using dot properties syntax' () {

        given:
        def folder = Files.createTempDirectory('test')
        def main = folder.resolve('nextflow.config')
        def snippet1 = folder.resolve('snippet1.config')
        def snippet2 = folder.resolve('snippet2.config')

        main.text = """
            process.name = 'alpha'

            includeConfig "$snippet1"
            """

        snippet1.text = """
            params.xxx = 'x'

            process.cpus = 4
            process.memory = '8GB'

            includeConfig "$snippet2"
            """

        snippet2.text = '''
            params.yyy = 'y'
            process { disk = '1TB' }
            process.resources.foo = 1
            process.resources.bar = 2
            '''

        when:
        def config = new ConfigParserV2().parse(main)
        then:
        config.params.xxx == 'x'
        config.params.yyy == 'y'
        config.process.name == 'alpha'
        config.process.cpus == 4
        config.process.memory == '8GB'
        config.process.disk == '1TB'
        config.process.resources.foo == 1
        config.process.resources.bar == 2

        cleanup:
        folder?.deleteDir()

    }

    def 'should parse multiple relative files' () {

        given:
        def folder = Files.createTempDirectory('test')
        def main = folder.resolve('main.config')
        def folder1 = folder.resolve('dir1')
        def folder2 = folder.resolve('dir2')
        def folder3 = folder.resolve('dir3')
        folder1.mkdirs()
        folder2.mkdirs()
        folder3.mkdirs()

        main. text = '''
            profiles {
                proc1 {
                    includeConfig 'dir1/config'
                }
                proc2 {
                    includeConfig 'dir2/config'
                }
                proc3 {
                    includeConfig 'dir3/config'
                }
            }
            '''

        folder.resolve('dir1/config').text = '''
            process.cpus = 4
            process.memory = '8GB'
            '''

        folder.resolve('dir2/config').text = '''
            process.cpus = 1
            process.memory = '6GB'
            '''

        folder.resolve('dir3/config').text = '''
            process.cpus = 32
            process.disk = '500GB'
            '''

        when:
        def config = new ConfigParserV2().parse(main)
        then:
        config.profiles.proc1.process.cpus == 4
        config.profiles.proc1.process.memory == '8GB'
        config.profiles.proc2.process.cpus == 1
        config.profiles.proc2.process.memory == '6GB'
        config.profiles.proc3.process.cpus == 32
        config.profiles.proc3.process.disk == '500GB'

        cleanup:
        folder?.deleteDir()

    }

    def 'should parse nested relative files' () {

        given:
        def folder = Files.createTempDirectory('test')
        def main = folder.resolve('main.config')
        def folder1 = folder.resolve('dir1/dir3')
        def folder2 = folder.resolve('dir2')
        folder1.mkdirs()
        folder2.mkdirs()

        main. text = """
            profiles {
                foo {
                    process.disk = '1TB'
                    includeConfig "dir1/nextflow.config"
                    includeConfig "dir2/nextflow.config"
                }
            }
            """

        folder.resolve('dir1/nextflow.config').text = """
            process.cpus = 4
            process.memory = '8GB'
            includeConfig 'dir3/nextflow.config'
            """

        folder.resolve('dir1/dir3/nextflow.config').text = '''
            process.ext.alpha = 1
            process.ext.delta = 2
            '''

        folder.resolve('dir2/nextflow.config').text = '''
            process.ext.cmd1 = 'echo true'
            process.ext.cmd2 = 'echo false'
            '''

        when:
        def config = new ConfigParserV2().parse(main)
        then:
        config.profiles.foo.process.cpus == 4
        config.profiles.foo.process.memory == '8GB'
        config.profiles.foo.process.disk == '1TB'

        config.profiles.foo.process.ext.alpha == 1
        config.profiles.foo.process.ext.delta == 2

        config.profiles.foo.process.ext.cmd1 == 'echo true'
        config.profiles.foo.process.ext.cmd2 == 'echo false'

        cleanup:
        folder?.deleteDir()

    }

    def 'should load selected profile configuration' () {

        given:
        def folder = Files.createTempDirectory('test')
        def main = folder.resolve('nextflow.config')
        def snippet1 = folder.resolve('snippet1.config')
        def snippet2 = folder.resolve('snippet2.config')

        snippet1.text = '''
            process {
                cpus =  1
                memory = '2GB'
                disk = '100GB'
            }
            '''

        snippet2.text = '''
            process {
                cpus =  8
                memory = '20GB'
                disk = '2TB'
            }
            '''

        main.text = """
            workDir = '/my/scratch'

            profiles {

                standard {}

                slow {
                    includeConfig "$snippet1"
                }

                fast {
                    workDir = '/fast/scratch'
                    includeConfig "$snippet2"
                }

            }
            """

        when:
        def config1 = new ConfigParserV2()
                        .setProfiles(['slow'])
                        .parse(main)
        then:
        config1.workDir == '/my/scratch'
        config1.process.cpus == 1
        config1.process.memory == '2GB'
        config1.process.disk == '100GB'

        when:
        def config2 = new ConfigParserV2()
                        .setProfiles(['fast'])
                        .parse(main)
        then:
        config2.workDir == '/fast/scratch'
        config2.process.cpus == 8
        config2.process.memory == '20GB'
        config2.process.disk == '2TB'

        cleanup:
        folder.deleteDir()

    }

    def 'should return the set of declared profiles' () {

        given:
        def text = '''
            profiles {
                alpha {
                    a = 1
                }
                beta {
                    b = 2
                }
            }
            '''

        when:
        def slurper = new ConfigParserV2().setProfiles(['alpha'])
        slurper.parse(text)
        then:
        slurper.getDeclaredProfiles() == ['alpha','beta'] as Set

        when:
        slurper = new ConfigParserV2().setProfiles(['omega'])
        slurper.parse(text)
        then:
        slurper.getDeclaredProfiles() == ['alpha','beta'] as Set
    }

    def 'should return the map of declared params' () {

        given:
        def text = '''
        params {
            a = 1
            b = 2
        }

        profiles {
            alpha {
                params.a = 3
            }
        }
        '''

        when:
        def slurper = new ConfigParserV2().setParams([c: 4])
        slurper.parse(text)
        then:
        slurper.getDeclaredParams() == [a: 1, b: 2]

        when:
        slurper = new ConfigParserV2().setParams([c: 4]).setProfiles(['alpha'])
        slurper.parse(text)
        then:
        slurper.getDeclaredParams() == [a: 3, b: 2]
    }

    def 'should ignore config includes when specified' () {
        given:
        def text = '''
            manifest {
                description = 'some text ..'
            }

            includeConfig 'this'
            includeConfig 'that'
            '''

        when:
        def config = new ConfigParserV2().setIgnoreIncludes(true).parse(text)
        then:
        config.manifest.description == 'some text ..'

        when:
        new ConfigParserV2().parse(text)
        then:
        thrown(NoSuchFileException)

    }

    def 'should parse file named as a top config scope' () {
        given:
        def folder = Files.createTempDirectory('test')
        def configFile = folder.resolve('XXX.config')
        configFile.text = 'XXX.enabled = true'

        when:
        new ConfigParserV2().parse(configFile)
        then:
        noExceptionThrown()

        cleanup:
        folder?.deleteDir()
    }

    def 'should access node metadata' () {

        given:
        def configText = '''
           params.str1 = 'hello'
           params.str2 = "${params.str1} world"
           process.clusterOptions = { "$params.str1" }
           process.ext = [foo: 'hello', bar: { params.str2 }]
           '''

        when:
        def config = new ConfigParserV2()
                .parse(configText)
        then:
        config.params.str1 instanceof String
        config.params.str2 instanceof GString
        config.process.clusterOptions instanceof Closure
        config.process.ext.bar instanceof Closure

        when:
        config = new ConfigParserV2()
                .setRenderClosureAsString(false)
                .parse(configText)
        then:
        config.params.str1 instanceof String
        config.params.str2 instanceof GString
        config.process.clusterOptions instanceof Closure
        config.process.ext.bar instanceof Closure

        when:
        config = new ConfigParserV2()
                .setRenderClosureAsString(true)
                .parse(configText)
        then:
        config.params.str1 == 'hello'
        config.params.str2 == 'hello world'
        config.process.clusterOptions == new ConfigClosurePlaceholder('{ "$params.str1" }')
        config.process.ext.foo == 'hello'
        config.process.ext.bar == new ConfigClosurePlaceholder('{ params.str2 }')

    }

    def 'should handle extend mem and duration units' () {
        ConfigObject result
        def CONFIG = '''
            mem1 = 1.GB
            mem2 = 1_000_000.toMemory()
            mem3 = MemoryUnit.of(2_000)
            time1 = 2.hours
            time2 = 60_000.toDuration()
            time3 = Duration.of(120_000)
            flag = 10000 < 1.GB
            '''

        when:
        result = new ConfigParserV2()
                .parse(CONFIG)
        then:
        result.mem1 instanceof MemoryUnit
        result.mem1 == MemoryUnit.of('1 GB')
        result.mem2 == MemoryUnit.of(1_000_000)
        result.mem3 == MemoryUnit.of(2_000)
        result.time1 instanceof Duration
        result.time1 == Duration.of('2 hours')
        result.time2 == Duration.of(60_000)
        result.time3 == Duration.of(120_000)
        result.flag == true
    }

    def 'should parse a config from an http server' () {
        given:
        def folder = Files.createTempDirectory('test')
        folder.resolve('conf').mkdir()

        HttpServer server = HttpServer.create(new InetSocketAddress(9900), 0);
        server.createContext("/", new ConfigFileHandler(folder));
        server.start()

        folder.resolve('nextflow.config').text = '''
            includeConfig 'conf/base.config'
            includeConfig 'http://localhost:9900/conf/remote.config'
            '''

        folder.resolve('conf/base.config').text = '''
            params.foo = 'Hello'
            params.bar = 'world!'
            '''

        folder.resolve('conf/remote.config').text = '''
            process {
                cpus = 4
                memory = '10GB'
            }
            '''

        when:
        def url = 'http://localhost:9900/nextflow.config' as Path
        def cfg = new ConfigBuilder().buildGivenFiles(url)
        then:
        cfg.params.foo == 'Hello'
        cfg.params.bar == 'world!'
        cfg.process.cpus == 4
        cfg.process.memory == '10GB'

        cleanup:
        server?.stop(0)
        folder?.deleteDir()
    }

    def 'should not overwrite param values with nested values of with the same name' () {
        given:
        def CONFIG = '''
            params {
              foo {
                bar = 'bar1'
              }
              baz = 'baz1'
              nested {
                baz = 'baz2'
                foo {
                  bar = 'bar2'
                }
              }
            }
            '''

        when:
        def config = new ConfigParserV2().parse(CONFIG)

        then:
        config.params.foo.bar == 'bar1'
        config.params.baz == 'baz1'
        config.params.nested.baz == 'baz2'
        config.params.nested.foo.bar == 'bar2'

        when:
        CONFIG = '''
            params {
              foo.bar = 'bar1'
              baz = 'baz1'
              nested.baz = 'baz2'
              nested.foo.bar = 'bar2'
            }
            '''
        and:
        config = new ConfigParserV2().parse(CONFIG)

        then:
        config.params.foo.bar == 'bar1'
        config.params.baz == 'baz1'
        config.params.nested.baz == 'baz2'
        config.params.nested.foo.bar == 'bar2'
    }

    def 'should apply profiles in the order they are specified at runtime' () {
        given:
        def CONFIG = '''
            profiles {
                foo {
                    params.input = 'foo'
                }

                bar {
                    params.input = 'bar'
                }
            }
            '''

        when:
        def config = new ConfigParserV2().setProfiles(['foo', 'bar']).parse(CONFIG)
        then:
        config.params.input == 'bar'

        when:
        config = new ConfigParserV2().setProfiles(['bar', 'foo']).parse(CONFIG)
        then:
        config.params.input == 'foo'
    }

    def 'should allow mixed use of dot and block syntax in a profile' () {
        given:
        def CONFIG = '''
            profiles {
                foo {
                    process.memory = '2 GB'
                    process {
                        cpus = 2
                    }
                }
            }
            '''

        when:
        def config = new ConfigParserV2().setProfiles(['foo']).parse(CONFIG)

        then:
        config.process.memory == '2 GB'
        config.process.cpus == 2
    }

    def 'should share params with included config files' () {

        given:
        def folder = Files.createTempDirectory('test')
        def main = folder.resolve('nextflow.config')
        def included = folder.resolve('included.config')

        main.text = """
            params.repo = 'foo/bar'

            includeConfig 'included.config'
            """

        included.text = '''
            params.url = "http://github.com/${params.repo}"
            '''

        when:
        def config = new ConfigParserV2().parse(main)
        then:
        config.params.repo == 'foo/bar'
        config.params.url == 'http://github.com/foo/bar'

        cleanup:
        folder?.deleteDir()

    }

    def 'should not evaluate profile blocks that were not selected' () {

        given:
        def folder = Files.createTempDirectory('test')
        def main = folder.resolve('nextflow.config')
        def included = folder.resolve('included.config')

        main.text = """
            profiles {
                bar {
                    includeConfig 'included.config'
                }
            }
            """

        included.text = '''
            // syntax error
            if( params )
                params.name = 'bar'
            '''

        when:
        new ConfigParserV2().setProfiles([]).parse(main)
        then:
        noExceptionThrown()

        when:
        new ConfigParserV2().setProfiles(['bar']).parse(main)
        then:
        thrown(ConfigParseException)

        when:
        new ConfigParserV2().parse(main)
        then:
        thrown(ConfigParseException)

        cleanup:
        folder?.deleteDir()

    }

    static class ConfigFileHandler implements HttpHandler {

        Path folder

        ConfigFileHandler(Path folder) {
             this.folder = folder
        }

        void handle(HttpExchange request) throws IOException {
            def path = request.requestURI.toString().substring(1)
            def file = folder.resolve(path)

            Headers header = request.getResponseHeaders()
            header.add("Content-Type", "text/plain")
            request.sendResponseHeaders(200, file.size())

            OutputStream os = request.getResponseBody();
            os.write(file.getBytes());
            os.close();
        }
    }

}
