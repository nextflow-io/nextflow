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

import nextflow.Global
import nextflow.Session
import nextflow.file.http.XPath
import nextflow.processor.TaskBean
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class FusionScriptLauncherTest extends Specification {

    def 'should get container mount' () {
        given:
        Global.session = Mock(Session) { getConfig() >> [:] }
        and:
        def fusion = new FusionScriptLauncher(scheme: 'http')

        when:
        def result = fusion.toContainerMount(XPath.get('http://foo/a/b/c.txt'))
        then:
        result == Path.of('/fusion/http/foo/a/b/c.txt')

        when:
        result = fusion.toContainerMount(XPath.get('http://foo/a/x/y.txt'))
        then:
        result == Path.of('/fusion/http/foo/a/x/y.txt')

        when:
        result = fusion.toContainerMount(XPath.get('http://bar/z.txt'))
        then:
        result == Path.of('/fusion/http/bar/z.txt')

    }

    def 'should get fusion env' () {
        given:
        Global.config = [:]
        and:
        def fusion = new FusionScriptLauncher(
                scheme: 'http',
                remoteWorkDir: XPath.get('http://foo/work'))

        expect:
        fusion.fusionEnv() == [
                FUSION_WORK: '/fusion/http/foo/work',
                FUSION_TAGS: "[.command.*|.exitcode|.fusion.*](nextflow.io/metadata=true),[*](nextflow.io/temporary=true)"
        ]
    }

    def 'should get fusion logs env' () {
        given:
        Global.config = [fusion: [logLevel:'debug', logOutput:'stdout', tags: false]]
        and:
        def fusion = new FusionScriptLauncher(
                scheme: 'http',
                remoteWorkDir: XPath.get('http://foo/work'))

        expect:
        fusion.fusionEnv() == [
                FUSION_WORK: '/fusion/http/foo/work',
                FUSION_LOG_LEVEL: 'debug',
                FUSION_LOG_OUTPUT: 'stdout'
        ]
    }

    def 'should get fusion with custom tags' () {
        given:
        Global.config = [fusion: [tags: 'custom-tags-pattern-here']]
        and:
        def fusion = new FusionScriptLauncher(
                scheme: 'http',
                remoteWorkDir: XPath.get('http://foo/work'))

        expect:
        fusion.fusionEnv() == [
                FUSION_WORK: '/fusion/http/foo/work',
                FUSION_TAGS: 'custom-tags-pattern-here'
        ]
    }

    def 'should get header script' () {
        given:
        def fusion = new FusionScriptLauncher(scheme: 's3')
        def task = Mock(TaskBean) { getWorkDir() >> Path.of('/some/work/dir')}

        expect:
        fusion.headerScript(task) == 'NXF_CHDIR=/some/work/dir\n'
    }
}
