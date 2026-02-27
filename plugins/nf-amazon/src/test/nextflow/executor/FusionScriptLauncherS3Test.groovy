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

package nextflow.executor

import java.nio.file.Path

import nextflow.Global
import nextflow.SysEnv
import nextflow.cloud.aws.util.S3PathFactory
import nextflow.fusion.FusionScriptLauncher
import nextflow.processor.TaskBean
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class FusionScriptLauncherS3Test extends Specification {

    def 'should get container mount' () {
        given:
        Global.config = Collections.emptyMap()
        and:
        def fusion = new FusionScriptLauncher(Mock(TaskBean), 's3', null)

        when:
        def result = fusion.toContainerMount(S3PathFactory.parse('s3://foo/a/b/c.txt'))
        then:
        result == Path.of('/fusion/s3/foo/a/b/c.txt')

        when:
        result = fusion.toContainerMount(S3PathFactory.parse('s3://foo/a/x/y.txt'))
        then:
        result == Path.of('/fusion/s3/foo/a/x/y.txt')

        when:
        result = fusion.toContainerMount(S3PathFactory.parse('s3://bar/z.txt'))
        then:
        result == Path.of('/fusion/s3/bar/z.txt')

    }


    def 'should get fusion env with s3 endpoint' () {
        given:
        Global.config = [:]
        and:
        SysEnv.push([AWS_S3_ENDPOINT: 'http://foo.com'])
        and:
        def fusion = new FusionScriptLauncher(Mock(TaskBean), 's3', S3PathFactory.parse('s3://foo/work'))

        expect:
        fusion.fusionEnv() == [AWS_S3_ENDPOINT: 'http://foo.com',
                               FUSION_WORK: '/fusion/s3/foo/work',
                               FUSION_TAGS: "[.command.*|.exitcode|.fusion.*](nextflow.io/metadata=true),[*](nextflow.io/temporary=true)"
                                ]

        cleanup:
        SysEnv.pop()
    }

    def 'should get fusion env with aws credentials' () {
        given:
        SysEnv.push([AWS_ACCESS_KEY_ID: 'xxx', AWS_SECRET_ACCESS_KEY: 'zzz'])
        and:
        Global.config = [fusion: [exportAwsAccessKeys: true]]
        and:
        def fusion = new FusionScriptLauncher(Mock(TaskBean), 's3', S3PathFactory.parse('s3://foo/work'))

        expect:
        fusion.fusionEnv() == [AWS_ACCESS_KEY_ID: 'xxx',
                               AWS_SECRET_ACCESS_KEY: 'zzz',
                               FUSION_WORK: '/fusion/s3/foo/work',
                               FUSION_TAGS: "[.command.*|.exitcode|.fusion.*](nextflow.io/metadata=true),[*](nextflow.io/temporary=true)"
        ]

        cleanup:
        Global.config = null
        SysEnv.pop()
    }

    def 'should get fusion env with aws credentials in nextflow config' () {
        given:
        SysEnv.push([:])
        and:
        def CONFIG = [fusion: [exportAwsAccessKeys: true], aws: [accessKey: 'k1', secretKey: 's1', client: [endpoint: 'http://minio.com']]]
        Global.config = CONFIG
        and:
        def fusion = new FusionScriptLauncher(Mock(TaskBean), 's3', S3PathFactory.parse('s3://foo/work'))

        expect:
        fusion.fusionEnv() == [AWS_ACCESS_KEY_ID: 'k1',
                               AWS_SECRET_ACCESS_KEY: 's1',
                               AWS_S3_ENDPOINT: 'http://minio.com',
                               FUSION_WORK: '/fusion/s3/foo/work',
                               FUSION_TAGS: "[.command.*|.exitcode|.fusion.*](nextflow.io/metadata=true),[*](nextflow.io/temporary=true)"
        ]

        cleanup:
        Global.config = null
        SysEnv.pop()
    }

}
