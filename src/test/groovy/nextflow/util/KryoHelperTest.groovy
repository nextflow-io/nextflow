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
import nextflow.container.ContainerConfig
import nextflow.file.FileHelper
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class KryoHelperTest extends  Specification {

    def testUUIDSerialization() {

        setup:
        def file = File.createTempFile('uuid-test',null)
        file.deleteOnExit()

        when:
        def uuid = UUID.randomUUID()
        KryoHelper.serialize(uuid, file)
        UUID copy = KryoHelper.deserialize(file)
        then:
        copy == uuid
    }


    def testUrlSerialization() {

        setup:
        def file = File.createTempFile('url-test',null)
        file.deleteOnExit()

        when:
        def url = file.toURI().toURL()
        KryoHelper.serialize(url, file)
        URL copy = KryoHelper.deserialize(file)
        then:
        copy == url
    }

    def testBytes() {

        setup:
        def x = 'Hola'

        when:
        def buffer = KryoHelper.serialize(x)
        def copy = KryoHelper.deserialize(buffer)

        then:
        copy == x

    }


    def testGString() {
        setup:
        def w = "world!"
        def x = "Hello $w"

        when:
        def buffer = KryoHelper.serialize(x)
        def copy = KryoHelper.deserialize(buffer)

        then:
        copy == x
    }

    def testDuration() {
        given:
        def d = Duration.of('24h')
        when:
        def buffer = KryoHelper.serialize(d)
        then:
        KryoHelper.deserialize(buffer) == d
        KryoHelper.deserialize(buffer).toString() == '1d'
    }

    def testMemUnit() {
        given:
        def m = new MemoryUnit('100MB')
        when:
        def buffer = KryoHelper.serialize(m)
        then:
        KryoHelper.deserialize(buffer) == m
        KryoHelper.deserialize(buffer).toString() == '100 MB'
    }

    def testSerializeArrayBag() {
        given:
        def bag = new ArrayBag(10,'Hello', 200F, 'World')
        when:
        def copy = KryoHelper.deserialize(KryoHelper.serialize(bag))
        then:
        copy == bag
        copy instanceof ArrayBag
        !copy.is(bag)
    }

    def testSerializeContainerConfig() {

        given:
        def cfg = new ContainerConfig([enabled: true, engine: 'docker', xxx: 'hello'])
        when:
        def copy = KryoHelper.deserialize(KryoHelper.serialize(cfg))
        then:
        copy == cfg
        copy instanceof ContainerConfig
        copy.engine == 'docker'
        copy.enabled == true
        copy.xxx == 'hello'

    }

    def 'should serialised a tuple array' () {

        given:
        def tuple = new ArrayTuple(['alpha','beta',null,'gamma'])
        when:
        def buffer = KryoHelper.serialize(tuple)
        then:
        KryoHelper.deserialize(buffer) == tuple
        KryoHelper.deserialize(buffer) instanceof ArrayTuple

    }

    def 'should serialise a map entry' () {

        given:
        def map = [foo:1]

        when:
        def entry = map.entrySet().iterator().next()
        def buffer = KryoHelper.serialize(entry)
        then:
        KryoHelper.deserialize(buffer) instanceof Map.Entry
        KryoHelper.deserialize(buffer) == entry

    }

    def 'should serialise xpath' () {
        when:
        def file = FileHelper.asPath('http://host.com/foo.txt')
        def buffer = KryoHelper.serialize(file)
        then:
        KryoHelper.deserialize(buffer).getClass().getName() == 'nextflow.file.http.XPath'
        KryoHelper.deserialize(buffer).toUri() == new URI('http://host.com/foo.txt')
    }


}
