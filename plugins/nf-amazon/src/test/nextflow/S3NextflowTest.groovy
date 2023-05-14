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

package nextflow

import java.nio.file.Paths

import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class S3NextflowTest extends Specification {


    def 'should return s3 uris'() {
        expect:
        Nextflow.file('s3://foo/data/file.log') == Paths.get(new URI('s3:///foo/data/file.log'))
    }


    def 'should resolve rel paths against env base' () {
        given:
        SysEnv.push(NXF_FILE_ROOT: 's3://some/base/dir')

        expect:
        Nextflow.file( 's3://abs/path/file.txt' ) == Paths.get(new URI('s3:///abs/path/file.txt'))
        and:
        Nextflow.file( 'file.txt' ) == Paths.get(new URI('s3:///some/base/dir/file.txt'))

        cleanup:
        SysEnv.pop()
    }

}
