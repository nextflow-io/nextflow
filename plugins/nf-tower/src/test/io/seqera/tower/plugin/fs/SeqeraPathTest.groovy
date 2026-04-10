/*
 * Copyright 2013-2026, Seqera Labs
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

package io.seqera.tower.plugin.fs

import io.seqera.tower.plugin.dataset.SeqeraDatasetClient
import spock.lang.Specification

/**
 * Unit tests for {@link SeqeraPath}.
 */
class SeqeraPathTest extends Specification {

    private SeqeraFileSystem mockFs() {
        def provider = new SeqeraFileSystemProvider()
        def client = Mock(SeqeraDatasetClient)
        return new SeqeraFileSystem(provider, client)
    }

    def "depth 0 - root path"() {
        given:
        def fs = mockFs()
        def path = new SeqeraPath(fs, 'seqera://')

        expect:
        path.depth() == 0
        path.isDirectory()
        !path.isRegularFile()
        path.org == null
        path.workspace == null
    }

    def "depth 1 - org path"() {
        given:
        def fs = mockFs()
        def path = new SeqeraPath(fs, 'seqera://acme')

        expect:
        path.depth() == 1
        path.isDirectory()
        !path.isRegularFile()
        path.org == 'acme'
        path.workspace == null
    }

    def "depth 2 - workspace path"() {
        given:
        def fs = mockFs()
        def path = new SeqeraPath(fs, 'seqera://acme/research')

        expect:
        path.depth() == 2
        path.isDirectory()
        path.org == 'acme'
        path.workspace == 'research'
        path.resourceType == null
    }

    def "depth 3 - resource type path"() {
        given:
        def fs = mockFs()
        def path = new SeqeraPath(fs, 'seqera://acme/research/datasets')

        expect:
        path.depth() == 3
        path.isDirectory()
        path.org == 'acme'
        path.workspace == 'research'
        path.resourceType == 'datasets'
        path.datasetName == null
    }

    def "depth 4 - dataset file path"() {
        given:
        def fs = mockFs()
        def path = new SeqeraPath(fs, 'seqera://acme/research/datasets/samples')

        expect:
        path.depth() == 4
        !path.isDirectory()
        path.isRegularFile()
        path.org == 'acme'
        path.workspace == 'research'
        path.resourceType == 'datasets'
        path.datasetName == 'samples'
        path.version == null
    }

    def "depth 4 - dataset with pinned version"() {
        given:
        def fs = mockFs()
        def path = new SeqeraPath(fs, 'seqera://acme/research/datasets/samples@2')

        expect:
        path.depth() == 4
        path.datasetName == 'samples'
        path.version == '2'
    }

    def "toUri round-trip - no version"() {
        given:
        def fs = mockFs()
        def uri = 'seqera://acme/research/datasets/samples'
        def path = new SeqeraPath(fs, uri)

        expect:
        path.toUri().toString() == uri
        path.toString() == uri
    }

    def "toUri round-trip - with version"() {
        given:
        def fs = mockFs()
        def uri = 'seqera://acme/research/datasets/samples@2'
        def path = new SeqeraPath(fs, uri)

        expect:
        path.toUri().toString() == uri
    }

    def "getParent - depth 4 returns depth 3"() {
        given:
        def fs = mockFs()
        def path = new SeqeraPath(fs, 'seqera://acme/research/datasets/samples')

        when:
        def parent = path.getParent()

        then:
        parent.toString() == 'seqera://acme/research/datasets'
        (parent as SeqeraPath).depth() == 3
    }

    def "getParent - depth 3 returns depth 2"() {
        given:
        def fs = mockFs()
        def path = new SeqeraPath(fs, 'seqera://acme/research/datasets')

        expect:
        path.getParent().toString() == 'seqera://acme/research'
    }

    def "getParent - depth 1 returns depth 0 root"() {
        given:
        def fs = mockFs()
        def path = new SeqeraPath(fs, 'seqera://acme')

        expect:
        path.getParent().toString() == 'seqera://'
    }

    def "getParent - root returns null"() {
        given:
        def fs = mockFs()
        def path = new SeqeraPath(fs, 'seqera://')

        expect:
        path.getParent() == null
    }

    def "resolve - appends segment to workspace"() {
        given:
        def fs = mockFs()
        def path = new SeqeraPath(fs, 'seqera://acme/research')

        when:
        def resolved = path.resolve('datasets')

        then:
        resolved.toString() == 'seqera://acme/research/datasets'
        (resolved as SeqeraPath).depth() == 3
    }

    def "resolve - appends dataset name to resource type"() {
        given:
        def fs = mockFs()
        def path = new SeqeraPath(fs, 'seqera://acme/research/datasets')

        when:
        def resolved = path.resolve('my-dataset')

        then:
        resolved.toString() == 'seqera://acme/research/datasets/my-dataset'
    }

    def "resolve - dataset name with version"() {
        given:
        def fs = mockFs()
        def path = new SeqeraPath(fs, 'seqera://acme/research/datasets')

        when:
        def resolved = path.resolve('samples@3')

        then:
        (resolved as SeqeraPath).datasetName == 'samples'
        (resolved as SeqeraPath).version == '3'
    }

    def "equality and hashCode"() {
        given:
        def fs = mockFs()
        def p1 = new SeqeraPath(fs, 'seqera://acme/research/datasets/samples')
        def p2 = new SeqeraPath(fs, 'seqera://acme/research/datasets/samples')
        def p3 = new SeqeraPath(fs, 'seqera://acme/research/datasets/other')

        expect:
        p1 == p2
        p1.hashCode() == p2.hashCode()
        p1 != p3
    }

    def "isAbsolute always true"() {
        given:
        def fs = mockFs()

        expect:
        new SeqeraPath(fs, 'seqera://acme').isAbsolute()
        new SeqeraPath(fs, 'seqera://').isAbsolute()
    }

    def "getNameCount equals depth"() {
        given:
        def fs = mockFs()

        expect:
        new SeqeraPath(fs, 'seqera://').nameCount == 0
        new SeqeraPath(fs, 'seqera://acme').nameCount == 1
        new SeqeraPath(fs, 'seqera://acme/research/datasets/samples').nameCount == 4
    }

    // ---- relativize ----

    def "relativize returns correct relative path string"() {
        given:
        def fs = mockFs()

        expect:
        new SeqeraPath(fs, base).relativize(new SeqeraPath(fs, other)).toString() == expected

        where:
        base                                            | other                                                        | expected
        'seqera://acme'                                 | 'seqera://acme/research'                                     | 'research'
        'seqera://acme/research'                        | 'seqera://acme/research/datasets'                            | 'datasets'
        'seqera://acme/research'                        | 'seqera://acme/research/datasets/samples'                    | 'datasets/samples'
        'seqera://acme/research/datasets'               | 'seqera://acme/research/datasets/samples'                    | 'samples'
        'seqera://acme/research/datasets/samples'       | 'seqera://acme/research/datasets/samples'                    | ''
    }

    def "relativize result round-trips through resolve"() {
        given:
        def fs = mockFs()
        def base = new SeqeraPath(fs, 'seqera://acme/research')
        def target = new SeqeraPath(fs, 'seqera://acme/research/datasets/samples')

        when:
        def rel = base.relativize(target)
        def restored = base.resolve(rel)

        then:
        rel.toString() == 'datasets/samples'
        !rel.isAbsolute()
        restored == target
    }

    def "relativize throws IllegalArgumentException when other does not start with this"() {
        given:
        def fs = mockFs()

        when:
        new SeqeraPath(fs, 'seqera://acme/research').relativize(new SeqeraPath(fs, 'seqera://other/ws'))

        then:
        thrown(IllegalArgumentException)
    }

    // ---- multi-segment resolve ----

    def "resolve with multi-segment string builds correct path"() {
        given:
        def fs = mockFs()
        def base = new SeqeraPath(fs, 'seqera://acme/research')

        expect:
        base.resolve('datasets/samples').toString() == 'seqera://acme/research/datasets/samples'
        base.resolve('datasets').toString() == 'seqera://acme/research/datasets'
    }

    def "resolve with absolute seqera URI returns that URI"() {
        given:
        def fs = mockFs()
        def base = new SeqeraPath(fs, 'seqera://acme/research')
        def absolute = 'seqera://other/ws/datasets/report'

        expect:
        base.resolve(absolute).toString() == absolute
    }

    def "isAbsolute is false for relative paths produced by relativize"() {
        given:
        def fs = mockFs()
        def rel = new SeqeraPath(fs, 'seqera://acme').relativize(new SeqeraPath(fs, 'seqera://acme/research'))

        expect:
        !rel.isAbsolute()
        rel.toString() == 'research'
    }

    // ---- getFileName ----

    def "getFileName returns relative path for each depth"() {
        given:
        def fs = mockFs()

        expect:
        new SeqeraPath(fs, 'seqera://').getFileName() == null
        new SeqeraPath(fs, 'seqera://acme').getFileName().toString() == 'acme'
        !new SeqeraPath(fs, 'seqera://acme').getFileName().isAbsolute()
        new SeqeraPath(fs, 'seqera://acme/research').getFileName().toString() == 'research'
        new SeqeraPath(fs, 'seqera://acme/research/datasets').getFileName().toString() == 'datasets'
        new SeqeraPath(fs, 'seqera://acme/research/datasets/samples').getFileName().toString() == 'samples'
        new SeqeraPath(fs, 'seqera://acme/research/datasets/samples@2').getFileName().toString() == 'samples@2'
    }

    def "getFileName is not absolute (uses relative constructor)"() {
        given:
        def fs = mockFs()
        def name = new SeqeraPath(fs, 'seqera://acme/research').getFileName()

        expect:
        !name.isAbsolute()
        name.toString() == 'research'
        name.getFileSystem() == null
    }

    // ---- asUri ----

    def "asUri - valid full path round-trips"() {
        expect:
        SeqeraPath.asUri('seqera://acme/research/datasets/samples').toString() == 'seqera://acme/research/datasets/samples'
        SeqeraPath.asUri('seqera://acme/research').toString() == 'seqera://acme/research'
    }

    def "asUri - empty path returns root URI"() {
        expect:
        SeqeraPath.asUri('seqera://').toString() == 'seqera:///'
    }

    def "asUri - path starting with dot has dot stripped"() {
        expect:
        // seqera://. → strips dot → seqera:// → hits empty-path case → seqera:///
        SeqeraPath.asUri('seqera://.').toString() == 'seqera:///'
        // seqera://./foo/bar → strips dot only (substring from index 10) → seqera:///foo/bar
        SeqeraPath.asUri('seqera://./foo/bar').toString() == 'seqera://foo/bar'
    }

    def "asUri - triple slash path throws IllegalArgumentException"() {
        when:
        SeqeraPath.asUri('seqera:///something')

        then:
        thrown(IllegalArgumentException)
    }

    def "asUri - missing protocol prefix throws IllegalArgumentException"() {
        when:
        SeqeraPath.asUri('s3://bucket/key')

        then:
        thrown(IllegalArgumentException)
    }

    def "asUri - null or empty throws IllegalArgumentException"() {
        when:
        SeqeraPath.asUri(null)

        then:
        thrown(IllegalArgumentException)

        when:
        SeqeraPath.asUri('')

        then:
        thrown(IllegalArgumentException)
    }

    // ---- startsWith ----

    def "startsWith - same path returns true"() {
        given:
        def fs = mockFs()
        def path = new SeqeraPath(fs, 'seqera://acme/research/datasets/samples')

        expect:
        path.startsWith(new SeqeraPath(fs, 'seqera://acme/research/datasets/samples'))
    }

    def "startsWith - prefix path returns true"() {
        given:
        def fs = mockFs()
        def path = new SeqeraPath(fs, 'seqera://acme/research/datasets/samples')

        expect:
        path.startsWith(new SeqeraPath(fs, 'seqera://acme'))
        path.startsWith(new SeqeraPath(fs, 'seqera://acme/research'))
        path.startsWith(new SeqeraPath(fs, 'seqera://'))
    }

    def "startsWith - component-wise not substring"() {
        given:
        def fs = mockFs()
        def path = new SeqeraPath(fs, 'seqera://acme-corp/research/datasets/samples')

        expect: 'acme is a substring prefix of acme-corp but not a component prefix'
        !path.startsWith(new SeqeraPath(fs, 'seqera://acme'))
    }

    def "startsWith - longer path returns false"() {
        given:
        def fs = mockFs()
        def path = new SeqeraPath(fs, 'seqera://acme')

        expect:
        !path.startsWith(new SeqeraPath(fs, 'seqera://acme/research'))
    }

    def "startsWith with string"() {
        given:
        def fs = mockFs()
        def path = new SeqeraPath(fs, 'seqera://acme/research/datasets/samples')

        expect:
        path.startsWith('seqera://acme')
        !path.startsWith('seqera://acm')
    }

    // ---- endsWith ----

    def "endsWith - absolute path requires exact match"() {
        given:
        def fs = mockFs()
        def path = new SeqeraPath(fs, 'seqera://acme/research/datasets/samples')

        expect:
        path.endsWith(new SeqeraPath(fs, 'seqera://acme/research/datasets/samples'))
        !path.endsWith(new SeqeraPath(fs, 'seqera://acme/research'))
    }

    def "endsWith - relative path matches trailing components"() {
        given:
        def fs = mockFs()
        def path = new SeqeraPath(fs, 'seqera://acme/research/datasets/samples')

        expect:
        path.endsWith(new SeqeraPath('samples'))
        path.endsWith(new SeqeraPath('datasets/samples'))
        !path.endsWith(new SeqeraPath('other'))
    }

    def "endsWith - component-wise not substring"() {
        given:
        def fs = mockFs()
        def path = new SeqeraPath(fs, 'seqera://acme/research/datasets/my-samples')

        expect: 'samples is a substring suffix of my-samples but not a component match'
        !path.endsWith(new SeqeraPath('samples'))
    }

    // ---- iterator ----

    def "iterator returns relative name components"() {
        given:
        def fs = mockFs()
        def path = new SeqeraPath(fs, 'seqera://acme/research/datasets/samples')

        when:
        def parts = path.iterator().collect { it.toString() }

        then:
        parts == ['acme', 'research', 'datasets', 'samples']
        path.iterator().every { !it.isAbsolute() }
    }

    def "iterator on root returns empty"() {
        given:
        def fs = mockFs()
        def path = new SeqeraPath(fs, 'seqera://')

        expect:
        !path.iterator().hasNext()
    }

    def "iterator on org returns single element"() {
        given:
        def fs = mockFs()
        def path = new SeqeraPath(fs, 'seqera://acme')

        when:
        def parts = path.iterator().collect { it.toString() }

        then:
        parts == ['acme']
    }
}
