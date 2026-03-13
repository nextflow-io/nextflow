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

package io.seqera.tower.plugin.dataset

import java.nio.file.AccessMode
import java.nio.file.FileSystemAlreadyExistsException
import java.nio.file.ReadOnlyFileSystemException

import spock.lang.Unroll
import spock.lang.Specification

/**
 * @author Edmund Miller
 */
class DatasetFileSystemProviderTest extends Specification {

    def 'should return dataset scheme'() {
        expect:
        new DatasetFileSystemProvider().getScheme() == 'dataset'
    }

    def 'should create filesystem'() {
        given:
        def provider = new DatasetFileSystemProvider()

        when:
        def fs = provider.newFileSystem(new URI('dataset:///'), [:])

        then:
        fs instanceof DatasetFileSystem
        fs.isOpen()
        fs.isReadOnly()
    }

    def 'should throw on duplicate filesystem creation'() {
        given:
        def provider = new DatasetFileSystemProvider()
        provider.newFileSystem(new URI('dataset:///'), [:])

        when:
        provider.newFileSystem(new URI('dataset:///'), [:])

        then:
        thrown(FileSystemAlreadyExistsException)
    }

    def 'should get path from URI'() {
        given:
        def provider = new DatasetFileSystemProvider()

        when:
        def path = provider.getPath(new URI('dataset://my-data'))

        then:
        path instanceof DatasetPath
        (path as DatasetPath).datasetName == 'my-data'
    }

    @Unroll
    def 'should reject non-dataset URI #uriString'() {
        given:
        def provider = new DatasetFileSystemProvider()

        when:
        provider.getPath(new URI(uriString))

        then:
        thrown(IllegalArgumentException)

        where:
        uriString << ['s3://bucket/key', 'gs://bucket/key', 'file:///tmp/data.csv']
    }

    @Unroll
    def 'should throw ReadOnlyFileSystemException on #operation'() {
        given:
        def provider = new DatasetFileSystemProvider()

        when:
        invoke.call(provider)

        then:
        thrown(ReadOnlyFileSystemException)

        where:
        operation             | invoke
        'createDirectory'     | { DatasetFileSystemProvider p -> p.createDirectory(p.getPath(new URI('dataset://test'))) }
        'delete'              | { DatasetFileSystemProvider p -> p.delete(p.getPath(new URI('dataset://test'))) }
        'copy'                | { DatasetFileSystemProvider p -> p.copy(p.getPath(new URI('dataset://src')), p.getPath(new URI('dataset://dst'))) }
        'move'                | { DatasetFileSystemProvider p -> p.move(p.getPath(new URI('dataset://src')), p.getPath(new URI('dataset://dst'))) }
        'write access check'  | { DatasetFileSystemProvider p -> p.checkAccess(p.getPath(new URI('dataset://test')), AccessMode.WRITE) }
    }

    def 'should not be hidden'() {
        given:
        def provider = new DatasetFileSystemProvider()
        def path = provider.getPath(new URI('dataset://test'))

        expect:
        !provider.isHidden(path)
    }

    @Unroll
    def 'should report isSameFile=#expected for #left vs #right'() {
        given:
        def provider = new DatasetFileSystemProvider()
        def a = provider.getPath(new URI(left))
        def b = provider.getPath(new URI(right))

        expect:
        provider.isSameFile(a, b) == expected

        where:
        left                | right               | expected
        'dataset://test'    | 'dataset://test'    | true
        'dataset://test1'   | 'dataset://test2'   | false
    }
}
