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

package nextflow.fusion

import java.nio.file.Path

import nextflow.container.ContainerConfig
import nextflow.file.http.XPath
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class FusionHelperTest extends Specification {

    def 'should make foreign path to fusion paths' () {
        when:
        def result = FusionHelper.toContainerMount(XPath.get('http://foo/a/b/c.txt'), 'http')
        then:
        result == Path.of('/fusion/http/foo/a/b/c.txt')

        when:
        result = FusionHelper.toContainerMount(XPath.get('http://foo/a/x/y.txt'), 'http')
        then:
        result == Path.of('/fusion/http/foo/a/x/y.txt')

        when:
        result = FusionHelper.toContainerMount(XPath.get('http://bar/z.txt'), 'http')
        then:
        result == Path.of('/fusion/http/bar/z.txt')

    }

    def 'should return fusion container command' () {
        given:
        def launcher = Mock(FusionScriptLauncher) {
            getEnvironment() >> ENV
        }

        when:
        def result = FusionHelper.runWithContainer(launcher, new ContainerConfig(CONFIG), NAME, OPTS, CMD)
        then:
        1 * launcher.fusionEnv() >> ENV
        and:
        result == EXPECTED

        where:
        CONFIG                  | ENV               | NAME          | OPTS          | CMD                   | EXPECTED
        [engine:'docker']       | [:]               | 'image:1'     | null          | ['echo', 'hello']     | "docker run -i --rm --privileged image:1 echo 'hello'"
        [engine:'docker']       | [FOO:'one']       | 'image:2'     | null          | ['echo', 'hello']     | "docker run -i -e \"FOO=one\" --rm --privileged image:2 echo 'hello'"
        [engine:'docker']       | [FOO:'one']       | 'image:2'     | '--this=that' | ['echo', 'hello']     | "docker run -i -e \"FOO=one\" --this=that --rm --privileged image:2 echo 'hello'"
        and:
        [engine:'singularity']  | [:]               | 'image:1'     | null          | ['echo', 'hello']     | "set +u; env - PATH=\"\$PATH\" \${TMP:+SINGULARITYENV_TMP=\"\$TMP\"} \${TMPDIR:+SINGULARITYENV_TMPDIR=\"\$TMPDIR\"} singularity exec --no-home --pid image:1 echo 'hello'"
        [engine:'singularity']  | [FOO:'one']       | 'image:1'     | null          | ['echo', 'hello']     | "set +u; env - PATH=\"\$PATH\" \${TMP:+SINGULARITYENV_TMP=\"\$TMP\"} \${TMPDIR:+SINGULARITYENV_TMPDIR=\"\$TMPDIR\"} SINGULARITYENV_FOO=\"one\" singularity exec --no-home --pid image:1 echo 'hello'"
        [engine:'singularity']  | [FOO:'one']       | 'image:1'     | '--this=that' | ['echo', 'hello']     | "set +u; env - PATH=\"\$PATH\" \${TMP:+SINGULARITYENV_TMP=\"\$TMP\"} \${TMPDIR:+SINGULARITYENV_TMPDIR=\"\$TMPDIR\"} SINGULARITYENV_FOO=\"one\" singularity exec --no-home --pid --this=that image:1 echo 'hello'"

    }

}
