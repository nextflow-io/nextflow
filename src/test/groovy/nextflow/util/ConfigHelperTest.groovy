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

package nextflow.util

import java.nio.file.Files
import java.nio.file.Paths

import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ConfigHelperTest extends Specification {

    @Unroll
    def "should parse string value: #str" () {

        expect:
        ConfigHelper.parseValue(str) == value

        where:
        str         | value
        'hola'      | 'hola'
        '1'         | 1
        "${Long.MAX_VALUE}" | Long.MAX_VALUE
        'True'      | true
        'False'     | false
        "10.2"      | 10.2
        '5sec'      | Duration.of('5sec')
        'live_in_3d'| 'live_in_3d'

    }

    def testResolveClasspaths() {

        given:
        def path1 = Files.createTempDirectory('path1')
        path1.resolve('file1').text = 'File 1'
        path1.resolve('file2.jar').text = 'File 2'
        path1.resolve('dir').mkdir()
        path1.resolve('dir/file3').text = 'File 3'
        path1.resolve('dir/file4').text = 'File 4'

        def path2 = Files.createTempDirectory('path2')
        path2.resolve('file5').text = 'File 5'
        path2.resolve('file6.jar').text = 'File 6'

        def path3 = Paths.get('/some/file')

        when:
        def list = ConfigHelper.resolveClassPaths([path1, path2, path3])
        then:
        list.size() == 4
        list.contains( path1 )
        list.contains( path1.resolve('file2.jar') )
        list.contains( path2 )
        list.contains( path2.resolve('file6.jar') )

        cleanup:
        path1?.deleteDir()
        path2?.deleteDir()

    }

    def 'should properties notation' () {

        given:
        def mem = { return 0 }
        def config = new ConfigObject()
        config.docker.enabled = true
        config.process.executor = 'slurm'
        config.process.queue = 'long'
        config.process.memory = mem
        config.process.omega = "Hi' there"

        when:
        def result = ConfigHelper.toPropertiesString(config)

        then:
        result == """
                docker.enabled=true
                process.omega=Hi' there
                process.queue=long
                process.memory=${mem.toString()}
                process.executor=slurm
                """
                .stripIndent().leftTrim()

    }
    def 'should canonical notation' () {

        given:
        def config = new ConfigObject()
        config.process.executor = 'slurm'
        config.process.queue = 'long'
        config.docker.enabled = true
        config.dummy = new ConfigObject() // <-- empty config object should not be print
        config.mail.from = 'yo@mail.com'
        config.mail.smtp.port = 25
        config.mail.smtp.user = 'yo'
        config.mail.smtp.host = 'mail.com'

        when:
        def result = ConfigHelper.toCanonicalString(config)
        then:
        result == '''
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
                          port = 25
                          user = 'yo'
                          host = 'mail.com'
                       }
                    }
                    '''
                .stripIndent().leftTrim()


        when:
        result = ConfigHelper.toCanonicalString(config, true)
        then:
        result == '''
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

    def 'should quote special identifiers' () {

        given:
        def config = new ConfigObject()
        config.'alpha-bravo' = 1
        config.'fizz+buzz'.x = 'hello'
        config.'fizz+buzz'.y = 'world'

        when:
        def result = ConfigHelper.toCanonicalString(config)
        then:
        result == '''
                    'alpha-bravo' = 1

                    'fizz+buzz' {
                       x = 'hello'
                       y = 'world'
                    }
                   '''
                    .stripIndent().leftTrim()
    }

    def 'should flatten notation' () {

        given:
        def config = new ConfigObject()
        config.process.queue = 'long'
        config.process.executor = 'slurm'
        config.docker.enabled = true
        config.zeta.'quoted-attribute'.foo = 1

        when:
        def result = ConfigHelper.toFlattenString(config)
        then:
        result == '''
                process.queue = 'long'
                process.executor = 'slurm'
                docker.enabled = true
                zeta.'quoted-attribute'.foo = 1
                '''
                .stripIndent().leftTrim()

        when:
        result = ConfigHelper.toFlattenString(config, true)
        then:
        result == '''
                docker.enabled = true
                process.executor = 'slurm'
                process.queue = 'long'
                zeta.'quoted-attribute'.foo = 1
                '''
                .stripIndent().leftTrim()

        when:
        config = new ConfigObject()
        config.foo = "Hi' there"

        result = ConfigHelper.toFlattenString(config)
        then:
        result == "foo = 'Hi\\' there'\n"

    }


    def 'should sort property keys' () {

        given:
        def props = new ConfigHelper.OrderedProperties()
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
        def props = new ConfigHelper.OrderedProperties(config)
        then:
        props.'omega' == 1
        props.'alpha' == 4
        props.'delta.y' == 'Hello'
        props.'delta.z' == 'world'
        props.keys().toSet() == ['alpha', 'delta.y', 'delta.z','omega'] as Set
    }

    def 'should verify valid identifiers' () {

        expect:
        ConfigHelper.isValidIdentifier('foo')
        ConfigHelper.isValidIdentifier('foo_bar')
        !ConfigHelper.isValidIdentifier('foo-bar')
        !ConfigHelper.isValidIdentifier('foo+bar')
        !ConfigHelper.isValidIdentifier('0foo')
    }

    @Unroll
    def 'should wrap names with key #str' () {
        expect:
        ConfigHelper.wrap0(str) == expected0
        ConfigHelper.wrap1(str) == expected1

        where:
        str                 | expected0             | expected1
        "foo"               | "foo"                 | 'foo'
        "foo-1"             | "'foo-1'"             | "'foo-1'"

        "withLabel:foo"     | "'withLabel:foo'"     | "withLabel:foo"
        "withLabel:1foo"    | "'withLabel:1foo'"    | "withLabel:'1foo'"

        "withName:foo"      | "'withName:foo'"      | "withName:foo"
        "withName:2foo"     | "'withName:2foo'"     | "withName:'2foo'"
    }



}
