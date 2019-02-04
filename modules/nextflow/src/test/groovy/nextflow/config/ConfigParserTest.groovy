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

package nextflow.config

import java.nio.file.Files
import java.nio.file.NoSuchFileException
import java.nio.file.Path

import com.sun.net.httpserver.Headers
import com.sun.net.httpserver.HttpExchange
import com.sun.net.httpserver.HttpHandler
import com.sun.net.httpserver.HttpServer
import spock.lang.Specification

import nextflow.util.Duration
import nextflow.util.MemoryUnit

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ConfigParserTest extends Specification {

    def 'should parse composed config files' () {

        given:
        def folder = File.createTempDir()
        def snippet1 = new File(folder,'config1.txt').absoluteFile
        def snippet2 = new File(folder,'config2.txt').absoluteFile


        def text = """
        process {
            name = 'alpha'
            resources {
                disk = '1TB'
                includeConfig "$snippet1"
            }
        }
        """

        snippet1.text = """
        cpus = 4
        memory = '8GB'
        nested {
            includeConfig("$snippet2")
        }
        """

        snippet2.text = '''
        foo = 1
        bar = 2
        '''

        when:
        def config = new ConfigParser().parse(text)
        then:
        config.process.name == 'alpha'
        config.process.resources.cpus == 4
        config.process.resources.memory == '8GB'
        config.process.resources.disk == '1TB'
        config.process.resources.nested.foo == 1
        config.process.resources.nested.bar == 2

        when:
        def buffer = new StringWriter()
        config.writeTo(buffer)
        def str = buffer.toString()
        then:
        str == '''
        process {
        	name='alpha'
        	resources {
        		disk='1TB'
        		cpus=4
        		memory='8GB'
        		nested {
        			foo=1
        			bar=2
        		}
        	}
        }
        '''.stripIndent().leftTrim()

        cleanup:
        folder?.deleteDir()

    }

    def 'should parse include config using dot properties syntax' () {

        given:
        def folder = File.createTempDir()
        def snippet1 = new File(folder,'config1.txt').absoluteFile
        def snippet2 = new File(folder,'config2.txt').absoluteFile


        def text = """
        process.name = 'alpha'
        includeConfig "$snippet1"
        """

        snippet1.text = """
        params.xxx = 'x'

        process.cpus = 4
        process.memory = '8GB'

        includeConfig("$snippet2")
        """

        snippet2.text = '''
        params.yyy = 'y'
        process { disk = '1TB' }
        process.resources.foo = 1
        process.resources.bar = 2
        '''

        when:
        def config = new ConfigParser().setBinding().parse(text)
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
        def folder = File.createTempDir()
        def main = new File(folder, 'main.config')
        def folder1 = new File(folder, 'dir1')
        def folder2 = new File(folder, 'dir2')
        def folder3 = new File(folder, 'dir3')
        folder1.mkdirs()
        folder2.mkdirs()
        folder3.mkdirs()

        main. text = '''
        profiles {
            includeConfig 'dir1/config'
            includeConfig 'dir2/config'
            includeConfig 'dir3/config'
        }
        '''

        new File(folder,'dir1/config').text = '''
        proc1 {
            cpus = 4
            memory = '8GB'
        }
        '''

        new File(folder, 'dir2/config').text = '''
        proc2 {
            cpus = MIN
            memory = '6GB'
        }
        '''

        new File(folder, 'dir3/config').text = '''
        proc3 {
            cpus = MAX
            disk = '500GB'
        }
        '''

        when:
        def config = new ConfigParser().setBinding([MIN: 1, MAX: 32]).parse(main)
        then:
        config.profiles.proc1.cpus == 4
        config.profiles.proc1.memory == '8GB'
        config.profiles.proc2.cpus == 1
        config.profiles.proc2.memory == '6GB'
        config.profiles.proc3.cpus == 32
        config.profiles.proc3.disk == '500GB'

        cleanup:
        folder?.deleteDir()

    }

    def 'should parse nested relative files' () {

        given:
        def folder = File.createTempDir()
        def main = new File(folder, 'main.config')
        def folder1 = new File(folder, 'dir1/dir3')
        def folder2 = new File(folder, 'dir2')
        folder1.mkdirs()
        folder2.mkdirs()

        main. text = """
        process {
            name = 'foo'
            resources {
                disk = '1TB'
                includeConfig "dir1/nextflow.config"
            }

            commands {
                includeConfig "dir2/nextflow.config"
            }
        }
        """

        new File(folder,'dir1/nextflow.config').text = """
        cpus = 4
        memory = '8GB'
        nested {
            includeConfig 'dir3/nextflow.config'
        }
        """

        new File(folder, 'dir1/dir3/nextflow.config').text = '''
        alpha = 1
        delta = 2
        '''

        new File(folder, 'dir2/nextflow.config').text = '''
        cmd1 = 'echo true'
        cmd2 = 'echo false'
        '''

        when:
        def config = new ConfigParser().parse(main)
        then:
        config.process.name == 'foo'
        config.process.resources.cpus == 4
        config.process.resources.memory == '8GB'
        config.process.resources.disk == '1TB'

        config.process.resources.nested.alpha == 1
        config.process.resources.nested.delta == 2

        config.process.commands.cmd1 == 'echo true'
        config.process.commands.cmd2 == 'echo false'

        cleanup:
        folder?.deleteDir()

    }



    def 'should load selected profile configuration' () {

        given:
        def folder = File.createTempDir()
        def snippet1 = new File(folder,'config1.txt').absoluteFile
        def snippet2 = new File(folder,'config2.txt').absoluteFile

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


        def configText = """
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
        def config1 = new ConfigParser()
                        .registerConditionalBlock('profiles','slow')
                        .parse(configText)
        then:
        config1.workDir == '/my/scratch'
        config1.process.cpus == 1
        config1.process.memory == '2GB'
        config1.process.disk == '100GB'

        when:
        def config2 = new ConfigParser()
                        .registerConditionalBlock('profiles','fast')
                        .parse(configText)
        then:
        config2.workDir == '/fast/scratch'
        config2.process.cpus == 8
        config2.process.memory == '20GB'
        config2.process.disk == '2TB'


        cleanup:
        folder.deleteDir()

    }

    def 'should return the set of visited block names' () {

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

        servers {
            local {
                x = 1
            }
            test {
                y = 2
            }
            prod {
                z = 3
            }
        }
        '''

        when:
        def slurper = new ConfigParser().registerConditionalBlock('profiles','alpha')
        slurper.parse(text)
        then:
        slurper.getConditionalBlockNames() == ['alpha','beta'] as Set

        when:
        slurper = new ConfigParser().registerConditionalBlock('profiles','omega')
        slurper.parse(text)
        then:
        slurper.getConditionalBlockNames() == ['alpha','beta'] as Set

        when:
        slurper = new ConfigParser().registerConditionalBlock('servers','xxx')
        slurper.parse(text)
        then:
        slurper.getConditionalBlockNames() == ['local','test','prod'] as Set

        when:
        slurper = new ConfigParser().registerConditionalBlock('foo','bar')
        slurper.parse(text)
        then:
        slurper.getConditionalBlockNames() == [] as Set
    }

    def 'should disable includeConfig parsing' () {
        given:
        def text = '''
        manifest {
          description = 'some text ..'
        }

        includeConfig 'this'
        includeConfig 'that'
        '''

        when:
        def config = new ConfigParser().setIgnoreIncludes(true).parse(text)
        then:
        config.manifest.description == 'some text ..'

        when:
        new ConfigParser().parse(text)
        then:
        thrown(NoSuchFileException)

    }

    def 'should parse file named as a top config scope' () {
        given:
        def folder = File.createTempDir()
        def configFile = new File(folder, 'XXX.config')
        configFile.text = 'XXX.enabled = true'

        when:
        new ConfigParser().parse(configFile)
        then:
        noExceptionThrown()

        cleanup:
        folder?.deleteDir()
    }


    def 'should access node metadata' () {

        given:
        ConfigObject result
        def CONFIG = '''   
           str1 = 'hello'
           str2 = "${str1} world"
           closure1 = { "$str" }
           map1 = [foo: 'hello', bar: { world }]
           '''

        when:
        result = new ConfigParser()
                .parse(CONFIG)
        then:
        result.str1 instanceof String
        result.str2 instanceof GString
        result.closure1 instanceof Closure
        result.map1.bar instanceof Closure

        when:
        result = new ConfigParser()
                .setRenderClosureAsString(false)
                .parse(CONFIG)
        then:
        result.str1 instanceof String
        result.str2 instanceof GString
        result.closure1 instanceof Closure
        result.map1.bar instanceof Closure

        when:
        result = new ConfigParser()
                    .setRenderClosureAsString(true)
                    .parse(CONFIG)
        then:
        result.str1 == 'hello'
        result.str2 == 'hello world'
        result.closure1 instanceof ConfigClosurePlaceholder
        result.closure1 == new ConfigClosurePlaceholder('{ "$str" }')
        result.map1.foo == 'hello'
        result.map1.bar == new ConfigClosurePlaceholder('{ world }')

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
        result = new ConfigParser()
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
        // launch web server
        HttpServer server = HttpServer.create(new InetSocketAddress(9900), 0);
        server.createContext("/", new ConfigFileHandler(folder));
        server.start()

        // main `nextflow.config` file
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
