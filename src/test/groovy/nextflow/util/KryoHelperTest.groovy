/*
 * Copyright (c) 2012, the authors.
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

import spock.lang.Ignore
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

    @Ignore
    def testClosureSerialization() {

        setup:
        def x = 1
        def y = 2
        def f = { return x+y }

        when:
        def buffer = KryoHelper.serialize(f)
        def copy = (Closure)KryoHelper.deserialize(buffer)

        then:
        copy.call() == 3

    }


}
