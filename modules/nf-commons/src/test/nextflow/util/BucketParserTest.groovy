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
