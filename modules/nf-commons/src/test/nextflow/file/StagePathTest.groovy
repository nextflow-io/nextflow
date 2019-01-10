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

package nextflow.file

import java.nio.file.Files

import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class StagePathTest extends Specification {

    def testToString() {

        when:
        def real = Files.createTempDirectory('test')
        def path = new StagePath(real)

        then:
        path.toString() == real.getFileName().toString()

        cleanup:
        real.deleteDir()

    }


    def testTextProperty() {

        given:
        def folder = Files.createTempDirectory('test')
        def real = folder.resolve('name')
        def stage = new StagePath(real)

        when:
        real.text = 'Hello'
        then:
        stage.text == 'Hello'

        when:
        stage.text = 'Ciao!'
        then:
        real.text == 'Ciao!'

        cleanup:
        folder.deleteDir()

    }

    def testReadLines() {

        given:
        def folder = Files.createTempDirectory('test')
        def real = folder.resolve('name')
        def stage = new StagePath(real)

        when:
        real.text = 'a\nbb\nccc'
        then:
        stage.readLines() == ['a','bb','ccc']

        cleanup:
        folder.deleteDir()

    }

    def testCopyTo() {

        given:
        def folder = Files.createTempDirectory('test')
        def target = Files.createTempDirectory('test')
        def real = folder.resolve('name')
        def stage = new StagePath(real)

        when:
        real.text = 'a\nbb\nccc'
        stage.copyTo(target)
        then:
        target.resolve('name').text == 'a\nbb\nccc'

        cleanup:
        folder.deleteDir()
        target.deleteDir()

    }
}
