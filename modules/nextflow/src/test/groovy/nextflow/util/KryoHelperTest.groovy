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
