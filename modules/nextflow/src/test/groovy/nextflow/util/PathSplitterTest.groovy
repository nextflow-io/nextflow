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
 *
 */

package nextflow.util

import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class PathSplitterTest extends Specification {

    @Unroll
    def 'should split path #PATH' () {
        expect:
        PathSplitter.parse(PATH) == RESULT

        where:
        PATH            | RESULT
        '/'             | new PathSplitter('/')
        'foo'           | new PathSplitter('foo', null)
        'foo/'          | new PathSplitter('foo', null)
        '/foo'          | new PathSplitter('/foo', null)
        '/foo/'         | new PathSplitter('/foo', null)
        and:
        'foo/bar/baz'   | new PathSplitter('foo', ['bar','baz'])
        '/foo/bar/baz'  | new PathSplitter('/foo', ['bar','baz'])
        'foo/bar/baz/'  | new PathSplitter('foo', ['bar','baz'])
        '/foo/bar/baz/' | new PathSplitter('/foo', ['bar','baz'])
        '/foo/x/y/z'    | new PathSplitter('/foo', ['x','y','z'])
        and:
        'file:/foo'                     | new PathSplitter('file:/foo', null)
        'file:/foo/x/y/z'               | new PathSplitter('file:/foo', ['x','y','z'])
        'file:///foo'                   | new PathSplitter('file:///foo', null)
        'file:///foo/x/y/z'             | new PathSplitter('file:///foo', ['x','y','z'])
        and:
        's3://my-bucket'                | new PathSplitter('s3://my-bucket')
        's3://my-bucket/'               | new PathSplitter('s3://my-bucket/')
        's3://my-bucket/foo'            | new PathSplitter('s3://my-bucket/foo')
        's3://my-bucket/foo/bar'        | new PathSplitter('s3://my-bucket/foo', ['bar'])
        's3://my-bucket/foo/bar/baz/'   | new PathSplitter('s3://my-bucket/foo', ['bar','baz'])
        and:
        'dx://my-bucket'                | new PathSplitter('dx://my-bucket')
        'dx://my-bucket:'               | new PathSplitter('dx://my-bucket:')
        'dx://my-bucket:/'              | new PathSplitter('dx://my-bucket:/')
        'dx://my-bucket:/foo'           | new PathSplitter('dx://my-bucket:/foo')
        'dx://my-bucket:/foo/bar'       | new PathSplitter('dx://my-bucket:/foo', ['bar'])
        'dx://my-bucket:/foo/bar/baz/'  | new PathSplitter('dx://my-bucket:/foo', ['bar','baz'])
    }
}
