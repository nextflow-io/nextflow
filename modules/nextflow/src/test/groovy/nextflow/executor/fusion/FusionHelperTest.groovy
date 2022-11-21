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

package nextflow.executor.fusion

import java.nio.file.Path

import nextflow.file.http.XPath
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class FusionHelperTest extends Specification {

    def 'should make foreign path to fusion paths' () {
        given:
        def buckets = new HashSet()

        when:
        def result = FusionHelper.toContainerMount(XPath.get('http://foo/a/b/c.txt'), 'http', buckets)
        then:
        result == Path.of('/fusion/http/foo/a/b/c.txt')
        and:
        buckets == [ 'foo' ] as Set

        when:
        result = FusionHelper.toContainerMount(XPath.get('http://foo/a/x/y.txt'), 'http', buckets)
        then:
        result == Path.of('/fusion/http/foo/a/x/y.txt')
        and:
        buckets == [ 'foo' ] as Set

        when:
        result = FusionHelper.toContainerMount(XPath.get('http://bar/z.txt'), 'http', buckets)
        then:
        result == Path.of('/fusion/http/bar/z.txt')
        and:
        buckets == [ 'foo', 'bar' ] as Set

    }

}
