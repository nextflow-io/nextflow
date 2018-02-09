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

import java.nio.file.NoSuchFileException

import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ComposedConfigSlurperTest extends Specification {

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
        def config = new ComposedConfigSlurper().parse(text)
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
        def config = new ComposedConfigSlurper().setBinding().parse(text)
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
        def config = new ComposedConfigSlurper().setBinding([MIN: 1, MAX: 32]).parse(main)
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
        def config = new ComposedConfigSlurper().parse(main)
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
        def config1 = new ComposedConfigSlurper()
                        .registerConditionalBlock('profiles','slow')
                        .parse(configText)
        then:
        config1.workDir == '/my/scratch'
        config1.process.cpus == 1
        config1.process.memory == '2GB'
        config1.process.disk == '100GB'

        when:
        def config2 = new ComposedConfigSlurper()
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
        def slurper = new ComposedConfigSlurper().registerConditionalBlock('profiles','alpha')
        slurper.parse(text)
        then:
        slurper.getConditionalBlockNames() == ['alpha','beta'] as Set

        when:
        slurper = new ComposedConfigSlurper().registerConditionalBlock('profiles','omega')
        slurper.parse(text)
        then:
        slurper.getConditionalBlockNames() == ['alpha','beta'] as Set

        when:
        slurper = new ComposedConfigSlurper().registerConditionalBlock('servers','xxx')
        slurper.parse(text)
        then:
        slurper.getConditionalBlockNames() == ['local','test','prod'] as Set

        when:
        slurper = new ComposedConfigSlurper().registerConditionalBlock('foo','bar')
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
        def config = new ComposedConfigSlurper().setIgnoreIncludes(true).parse(text)
        then:
        config.manifest.description == 'some text ..'

        when:
        new ComposedConfigSlurper().parse(text)
        then:
        thrown(NoSuchFileException)

    }

    def 'should parse file named as a top config scope' () {
        given:
        def folder = File.createTempDir()
        def configFile = new File(folder, 'XXX.config')
        configFile.text = 'XXX.enabled = true'

        when:
        new ComposedConfigSlurper().parse(configFile)
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
        result = new ComposedConfigSlurper()
                .parse(CONFIG)
        then:
        result.str1 instanceof String
        result.str2 instanceof GString
        result.closure1 instanceof Closure
        result.map1.bar instanceof Closure

        when:
        result = new ComposedConfigSlurper()
                .setRenderClosureAsString(false)
                .parse(CONFIG)
        then:
        result.str1 instanceof String
        result.str2 instanceof GString
        result.closure1 instanceof Closure
        result.map1.bar instanceof Closure

        when:
        result = new ComposedConfigSlurper()
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


}
