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

package nextflow.extension

import java.nio.file.Path

import com.google.cloud.storage.contrib.nio.CloudStorageFileSystem
import nextflow.util.Escape
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class EscapeTest2 extends Specification {


    Path asPath(String bucket, String path) {
        CloudStorageFileSystem.forBucket(bucket).getPath(path)
    }


    def 'should escape gs path'() {
        expect:
        Escape.uriPath(PATH) == EXPECTED
        where:
        PATH                                | EXPECTED
        asPath('foo','/work')   | 'gs://foo/work'
        asPath('foo','/b a r')  | 'gs://foo/b\\ a\\ r'
        asPath('f_o o','/bar')  | 'gs://f_o\\ o/bar'
    }
}
