/*
 * Copyright 2020-2022, Seqera Labs
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
 *
 */

package nextflow.file

import java.nio.file.Paths

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class FilePathFactoryTest extends Specification {

    def 'should parse file uris to path' () {
        given:
        def factory = new FilePathFactory()

        expect:
        factory.parseUri('file:/foo/bar') == Paths.get('/foo/bar')
        factory.parseUri('file:///foo/bar') == Paths.get('/foo/bar')
    }

    def 'should convert path to uri' () {
        given:
        def factory = new FilePathFactory()
        
        expect:
        factory.toUriString(Paths.get('/foo/bar')) == 'file:///foo/bar'
        factory.toUriString(Paths.get(new URI('file:/foo/bar'))) == 'file:///foo/bar'
        factory.toUriString(Paths.get(new URI('file:///foo/bar'))) == 'file:///foo/bar'
    }
}
