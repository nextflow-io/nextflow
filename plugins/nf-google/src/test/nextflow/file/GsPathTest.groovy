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

package nextflow.file

import com.google.cloud.storage.contrib.nio.CloudStoragePath
import nextflow.Global
import nextflow.Session
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class GsPathTest extends Specification {

    def 'should check equals and hashcode' () {
        given:
        Global.session = Mock(Session) {
            getConfig() >> [google:[project:'foo', region:'x']]
        }
        and:

        def path1 = FileHelper.asPath('gs://foo/some/foo.txt')
        def path2 = FileHelper.asPath('gs://foo/some/foo.txt')
        def path3 = FileHelper.asPath('gs://foo/some/bar.txt')
        def path4 = FileHelper.asPath('gs://bar/some/foo.txt')

        expect:
        path1 instanceof CloudStoragePath
        path2 instanceof CloudStoragePath
        path3 instanceof CloudStoragePath
        path4 instanceof CloudStoragePath

        and:
        path1 == path2
        path1 != path3
        path3 != path4

        and:
        path1.hashCode() == path2.hashCode()
        path1.hashCode() != path3.hashCode()
        path3.hashCode() != path4.hashCode()
    }

}
