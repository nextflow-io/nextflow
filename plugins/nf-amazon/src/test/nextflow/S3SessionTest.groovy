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
import nextflow.file.FileSystemPathFactory
import spock.lang.Unroll
import test.AppSpec
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class S3SessionTest extends AppSpec {

    private Path asPath(String str) {
        if( str==null )
            return null
        if( str.startsWith('s3://'))
            return FileSystemPathFactory.parse(str)
        else
            return Path.of(str)
    }

    @Unroll
    def 'should get cloud cache path' () {
        given:
        def session = Spy(Session)

        expect:
        session.cloudCachePath(CONFIG, FileHelper.asPath(WORKDIR)) == asPath(EXPECTED)

        where:
        CONFIG                                  | WORKDIR           | EXPECTED
        null                                    | '/foo'            | null
        [enabled:true]                          | 's3://foo/work'   | 's3://foo/work'
        [enabled:true, path:'s3://this/that']   | '/foo'            | 's3://this/that'

    }


    def 'should error with non-cloud bucket' () {
        given:
        def session = Spy(Session)

        when:
        session.cloudCachePath([enabled:true], Path.of('/foo/dir'))
        then:
        def e = thrown(IllegalArgumentException)
        e.message == "Storage path not supported by Cloud-cache - offending value: '/foo/dir'"

    }

}
