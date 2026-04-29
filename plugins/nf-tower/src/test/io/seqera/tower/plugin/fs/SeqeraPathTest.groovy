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

import spock.lang.Specification

/**
 * Unit tests for {@link SeqeraPath}.
 */
class SeqeraPathTest extends Specification {

    private SeqeraFileSystem mockFs() {
        def provider = new SeqeraFileSystemProvider()
        return new SeqeraFileSystem(provider)
    }

    // ---- depth / segment accessors ----

    def "depth 0 - root path"() {
        given:
        def fs = mockFs()
        def path = new SeqeraPath(fs, 'seqera://')

        expect:
        path.depth() == 0
        path.org == null
        path.workspace == null
        path.resourceType == null
        path.trail == []
    }

    def "depth 1 - org path"() {
        given:
        def fs = mockFs()
        def path = new SeqeraPath(fs, 'seqera://acme')

        expect:
        path.depth() == 1
        path.org == 'acme'
        path.workspace == null
    }

    def "depth 2 - workspace path"() {
        given:
        def fs = mockFs()
        def path = new SeqeraPath(fs, 'seqera://acme/research')

        expect:
        path.depth() == 2
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
        path.org == 'acme'
        path.workspace == 'research'
        path.resourceType == 'datasets'
        path.trail == []
    }

    def "depth 4 - dataset trail segment is raw (handler parses @version)"() {
        given:
        def fs = mockFs()
        def path = new SeqeraPath(fs, 'seqera://acme/research/datasets/samples')

        expect:
        path.depth() == 4
        path.resourceType == 'datasets'
        path.trail == ['samples']
    }

    def "dataset with @version suffix stays raw in trail"() {
        given:
        def fs = mockFs()
        def path = new SeqeraPath(fs, 'seqera://acme/research/datasets/samples@2')

        expect:
        // Path is resource-type-agnostic — no @version parsing here.
        path.depth() == 4
        path.trail == ['samples@2']
    }

    def "data-link path with provider, name, and sub-path"() {
        given:
        def fs = mockFs()
        def path = new SeqeraPath(fs, 'seqera://acme/research/data-links/AWS/inputs/reads/sample.fq.gz')

        expect:
        path.depth() == 7
        path.resourceType == 'data-links'
        path.trail == ['AWS', 'inputs', 'reads', 'sample.fq.gz']
    }

    def "data-link path at provider level"() {
        given:
        def fs = mockFs()
        def path = new SeqeraPath(fs, 'seqera://acme/research/data-links/AWS')

        expect:
        path.depth() == 4
        path.resourceType == 'data-links'
        path.trail == ['AWS']
    }

    // ---- toUri / toString ----

    def "toUri round-trip - no version"() {
        given:
        def fs = mockFs()
        def uri = 'seqera://acme/research/datasets/samples'
        def path = new SeqeraPath(fs, uri)

        expect:
        path.toUri().toString() == uri
        path.toString() == uri
    }

    def "toUri round-trip - dataset with @version"() {
        given:
        def fs = mockFs()
        def uri = 'seqera://acme/research/datasets/samples@2'
        def path = new SeqeraPath(fs, uri)

        expect:
        path.toUri().toString() == uri
    }

    def "toUri round-trip - deep data-link path"() {
        given:
        def fs = mockFs()
        def uri = 'seqera://acme/research/data-links/AWS/inputs/reads/sample.fq.gz'
        def path = new SeqeraPath(fs, uri)

        expect:
        path.toUri().toString() == uri
    }

    // ---- getParent ----

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

    def "getParent - depth 7 returns depth 6 (drops one trail segment)"() {
        given:
        def fs = mockFs()
        def path = new SeqeraPath(fs, 'seqera://acme/research/data-links/AWS/inputs/reads/s.fq.gz')

        when:
        def parent = path.getParent() as SeqeraPath

        then:
        parent.trail == ['AWS', 'inputs', 'reads']
        parent.depth() == 6
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

    // ---- resolve ----

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

    def "resolve - dataset name with @version preserved as raw trail segment"() {
        given:
        def fs = mockFs()
        def path = new SeqeraPath(fs, 'seqera://acme/research/datasets')

        when:
        def resolved = path.resolve('samples@3') as SeqeraPath

        then:
        resolved.trail == ['samples@3']
    }

    def "resolve - appends nested data-link path segment"() {
        given:
        def fs = mockFs()
        def path = new SeqeraPath(fs, 'seqera://acme/research/data-links/AWS/inputs')

        when:
        def child = path.resolve('reads') as SeqeraPath

        then:
        child.trail == ['AWS', 'inputs', 'reads']
    }

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

    // ---- equality / hashCode ----

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

    def "isAbsolute true when fs attached"() {
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
        new SeqeraPath(fs, 'seqera://acme/research/data-links/AWS/inputs/reads/a.fq').nameCount == 7
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

    def "relativize produces '..' segments for upward traversal"() {
        given:
        def fs = mockFs()

        expect:
        new SeqeraPath(fs, base).relativize(new SeqeraPath(fs, other)).toString() == expected

        where:
        base                                            | other                                            | expected
        'seqera://acme/research'                        | 'seqera://acme/dev'                              | '../dev'
        'seqera://acme/research/datasets'               | 'seqera://acme/dev'                               | '../../dev'
        'seqera://acme'                                 | 'seqera://other'                                  | '../other'
        'seqera://acme/ws1'                             | 'seqera://acme/ws2'                               | '../ws2'
        'seqera://acme/research/datasets/samples'       | 'seqera://acme/research/datasets/other'            | '../other'
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
        new SeqeraPath(fs, 'seqera://acme/research/data-links/AWS/inputs/reads/a.fq').getFileName().toString() == 'a.fq'
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
        SeqeraPath.asUri('seqera://.').toString() == 'seqera:///'
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

    // ---- trailing slash / accidental double-slash tolerance ----

    def "trailing slash on resource-type directory is ignored"() {
        given:
        def fs = mockFs()

        when:
        def p = new SeqeraPath(fs, 'seqera://acme/research/datasets/')

        then:
        p.depth() == 3
        p.resourceType == 'datasets'
        p.trail == []
    }

    def "trailing slash on data-link directory is ignored"() {
        given:
        def fs = mockFs()

        when:
        def p = new SeqeraPath(fs, 'seqera://acme/research/data-links/aws/inputs/')

        then:
        p.depth() == 5
        p.trail == ['aws', 'inputs']
    }

    def "accidental double-slash inside the trail is collapsed"() {
        given:
        def fs = mockFs()

        when:
        def p = new SeqeraPath(fs, 'seqera://acme/research/data-links/aws/inputs//reads/a.fq')

        then:
        p.trail == ['aws', 'inputs', 'reads', 'a.fq']
    }

    // ---- cached attributes ----

    def "cachedAttributes is null by default and preserved by resolveWithAttributes"() {
        given:
        def fs = mockFs()
        def parent = new SeqeraPath(fs, 'seqera://acme/research/data-links/aws/inputs')
        def attrs = new SeqeraFileAttributes(42L, java.time.Instant.EPOCH, java.time.Instant.EPOCH, 'k')

        when:
        def child = parent.resolveWithAttributes('reads', attrs)

        then:
        parent.cachedAttributes == null
        child.cachedAttributes === attrs
        child.toString() == 'seqera://acme/research/data-links/aws/inputs/reads'
    }

    def "cachedAttributes does not affect equals/hashCode"() {
        given:
        def fs = mockFs()
        def attrs = new SeqeraFileAttributes(true)
        def withAttrs = new SeqeraPath(fs, 'seqera://acme/research/datasets').resolveWithAttributes('samples', attrs)
        def withoutAttrs = new SeqeraPath(fs, 'seqera://acme/research/datasets/samples')

        expect:
        withAttrs == withoutAttrs
        withAttrs.hashCode() == withoutAttrs.hashCode()
    }

    def "iterator on deep data-link path returns all segments"() {
        given:
        def fs = mockFs()
        def path = new SeqeraPath(fs, 'seqera://acme/research/data-links/AWS/inputs/reads/a.fq')

        when:
        def parts = path.iterator().collect { it.toString() }

        then:
        parts == ['acme', 'research', 'data-links', 'AWS', 'inputs', 'reads', 'a.fq']
    }
}
