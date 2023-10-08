/*
 * Copyright 2013-2023, Seqera Labs
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

package nextflow

import java.nio.file.Path

import nextflow.file.FileHelper
import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class S3SessionTest extends Specification {

    @Unroll
    def 'should get cloud cache path' () {
        given:
        def session = Spy(Session)
        SysEnv.push(ENV)

        expect:
        session.initCloudCache(CONFIG, FileHelper.asPath(WORKDIR)) == EXPECTED

        cleanup:
        SysEnv.pop()

        where:
        CONFIG                                  | WORKDIR           | ENV           | EXPECTED
        null                                    | '/foo'            | [:]           | null
        [enabled:true, path:'s3://this/that']   | '/foo'            | [:]           | FileHelper.asPath('s3://this/that')
        [enabled:true, path:'s3://this/that']   | '/foo'            | [NXF_CLOUDCACHE_PATH:'s3://other/path']       | FileHelper.asPath('s3://this/that')
        [enabled:true]                          | '/foo'            | [NXF_CLOUDCACHE_PATH:'s3://other/path']       | FileHelper.asPath('s3://other/path')
        [enabled:true]                          | 's3://foo/work'   | [:]            | FileHelper.asPath('s3://foo/work')

    }


    def 'should error with non-cloud bucket' () {
        given:
        def session = Spy(Session)

        when:
        session.initCloudCache([enabled:true], Path.of('/foo/dir'))
        then:
        def e = thrown(IllegalArgumentException)
        e.message == "Storage path not supported by Cloud-cache - offending value: '/foo/dir'"

    }

}
