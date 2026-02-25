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

import spock.lang.Specification

/**
 * @author Edmund Miller
 */
class DatasetFileSystemProviderTest extends Specification {

    def 'should return dataset scheme'() {
        given:
        def provider = new DatasetFileSystemProvider()

        expect:
        provider.getScheme() == 'dataset'
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

    def 'should reject non-dataset URI'() {
        given:
        def provider = new DatasetFileSystemProvider()

        when:
        provider.getPath(new URI('s3://bucket/key'))

        then:
        thrown(IllegalArgumentException)
    }

    // -- read-only enforcement --

    def 'should throw ReadOnlyFileSystemException on createDirectory'() {
        given:
        def provider = new DatasetFileSystemProvider()
        def path = provider.getPath(new URI('dataset://test'))

        when:
        provider.createDirectory(path)

        then:
        thrown(ReadOnlyFileSystemException)
    }

    def 'should throw ReadOnlyFileSystemException on delete'() {
        given:
        def provider = new DatasetFileSystemProvider()
        def path = provider.getPath(new URI('dataset://test'))

        when:
        provider.delete(path)

        then:
        thrown(ReadOnlyFileSystemException)
    }

    def 'should throw ReadOnlyFileSystemException on copy'() {
        given:
        def provider = new DatasetFileSystemProvider()
        def src = provider.getPath(new URI('dataset://src'))
        def dst = provider.getPath(new URI('dataset://dst'))

        when:
        provider.copy(src, dst)

        then:
        thrown(ReadOnlyFileSystemException)
    }

    def 'should throw ReadOnlyFileSystemException on move'() {
        given:
        def provider = new DatasetFileSystemProvider()
        def src = provider.getPath(new URI('dataset://src'))
        def dst = provider.getPath(new URI('dataset://dst'))

        when:
        provider.move(src, dst)

        then:
        thrown(ReadOnlyFileSystemException)
    }

    def 'should throw ReadOnlyFileSystemException on write access check'() {
        given:
        def provider = new DatasetFileSystemProvider()
        def path = provider.getPath(new URI('dataset://test'))

        when:
        provider.checkAccess(path, AccessMode.WRITE)

        then:
        thrown(ReadOnlyFileSystemException)
    }

    def 'should not be hidden'() {
        given:
        def provider = new DatasetFileSystemProvider()
        def path = provider.getPath(new URI('dataset://test'))

        expect:
        !provider.isHidden(path)
    }

    def 'should detect same file for equal dataset paths'() {
        given:
        def provider = new DatasetFileSystemProvider()
        def a = provider.getPath(new URI('dataset://test'))
        def b = provider.getPath(new URI('dataset://test'))

        expect:
        provider.isSameFile(a, b)
    }

    def 'should detect different files for different dataset paths'() {
        given:
        def provider = new DatasetFileSystemProvider()
        def a = provider.getPath(new URI('dataset://test1'))
        def b = provider.getPath(new URI('dataset://test2'))

        expect:
        !provider.isSameFile(a, b)
    }
}
