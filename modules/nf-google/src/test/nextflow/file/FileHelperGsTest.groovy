/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

package nextflow.file

import java.nio.file.Paths

import com.google.cloud.storage.contrib.nio.CloudStorageFileSystem
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class FileHelperGsTest extends Specification {

    def 'should parse google storage path' () {

        expect:
        FileHelper.asPath('file.txt') ==
                Paths.get('file.txt')
        and:
        FileHelper.asPath('gs://foo') ==
                CloudStorageFileSystem.forBucket('foo').getPath('')

        and:
        FileHelper.asPath('gs://foo/this/and/that.txt') ==
                CloudStorageFileSystem.forBucket('foo').getPath('/this/and/that.txt')

        and:
        FileHelper.asPath('gs://foo/b a r.txt') ==
                CloudStorageFileSystem.forBucket('foo').getPath('/b a r.txt')
        
        and:
        FileHelper.asPath('gs://f o o/bar.txt') ==
                CloudStorageFileSystem.forBucket('f o o').getPath('/bar.txt')

        and:
        FileHelper.asPath('gs://f_o_o/bar.txt') ==
                CloudStorageFileSystem.forBucket('f_o_o').getPath('/bar.txt')
    }
}
