/*
 * Copyright (c) 2019-2020, Seqera Labs.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * This Source Code Form is "Incompatible With Secondary Licenses", as
 * defined by the Mozilla Public License, v. 2.0.
 */

package nextflow.util

import java.nio.file.Path

import nextflow.io.BucketParser
import spock.lang.Specification
import spock.lang.Unroll
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class BucketParserTest extends Specification {

    @Unroll
    def 'should parse bucket path #STR' () {
        expect:
        BucketParser.from(STR) == EXPECTED
        where:
        STR             | EXPECTED
        '/'             | new BucketParser(null, null, null)
        '/'             | new BucketParser(null, null, '/')
        '/a/b/c'        | new BucketParser(null, null, '/a/b/c')
        'a/b/c'         | new BucketParser(null, null, 'a/b/c')
        and:
        's3://foo'      | new BucketParser('s3','foo','/')
        's3://foo/'     | new BucketParser('s3','foo','/')
        's3://foo/x/y'  | new BucketParser('s3','foo','/x/y')
        's3://foo/x/y'  | new BucketParser('s3','foo','x/y')
        's3://foo/x/y/' | new BucketParser('s3','foo','/x/y')
        and:
        'gs://foo/x/y'  | new BucketParser('gs','foo','/x/y')
    }

    def 'should get key' () {
        expect:
        BucketParser.from(PATH).getKey() == KEY

        where:
        PATH                    | KEY
        's3://foo'              | ''
        's3://foo/bar/baz'      | 'bar/baz'
    }

    def 'should get tokens' () {
        given:
        def parsed = BucketParser.from(STR)
        expect:
        parsed.getBucket() == BUCKET
        parsed.getScheme() == SCHEME
        parsed.getPath() == Path.of(PATH)

        where:
        STR             | SCHEME    | BUCKET    | PATH
        '/'             | null      | null      | '/'
        '/a/b/c'        | null      | null      | '/a/b/c'
        'a/b/c'         | null      | null      | 'a/b/c'
        and:
        's3://foo'      | 's3'      | 'foo'     | '/'
        's3://foo'      | 's3'      | 'foo'     | '/'
        's3://foo/x/y'  | 's3'      | 'foo'     | '/x/y'
        's3://foo/x/y/' | 's3'      | 'foo'     | '/x/y'
        and:
        'gs://foo/x/y'  | 'gs'      | 'foo'     | '/x/y'

    }
}
