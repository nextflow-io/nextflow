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

package org.codehaus.groovy.runtime

import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths

import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class NioGroovyMethodsTest extends Specification {

    def testSize() {

        setup:
        Path path = Files.createTempFile('test_size', null);

        when:
        path.text = 'Hello'
        then:
        path.size() == 5

        cleanup:
        Files.delete(path)

    }

    static class Bean implements Serializable {
        def alpha
        def beta
        def delta
    }

    def testNewInputOutputStream() {

        setup:
        def bean = new Bean(alpha: 'a', beta: 2, delta: true)
        def path = Paths.get('file_data')

        when:
        def out = path.newObjectOutputStream()
        out.writeObject(bean)
        out.flush()
        def clone = (Bean)path.newObjectInputStream().readObject()
        then:
        clone.alpha == bean.alpha
        clone.beta == bean.beta
        clone.delta == bean.delta

        cleanup:
        path?.delete()

    }

}
