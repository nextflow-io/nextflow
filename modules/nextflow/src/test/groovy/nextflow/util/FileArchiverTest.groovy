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

package nextflow.util

import java.nio.file.Path

import nextflow.file.FileHelper
import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class FileArchiverTest extends Specification {

    def 'should empty target path' () {
        given:
        def helper = new FileArchiver([:])

        expect:
        helper.getBaseDir() == null
        helper.getTargetDir() == null
        and:
        helper.archivePath(null) == null
        helper.archivePath(Path.of('/some/data/file.txt')) == null
    }

    @Unroll
    def 'should final target path' () {
        given:
        def helper = new FileArchiver([NXF_ARCHIVE_DIR:'/data,http://bucket/data/export'])

        expect:
        helper.getBaseDir() == Path.of('/data')
        helper.getTargetDir() == FileHelper.asPath('http://bucket/data/export')

        and:
        helper.archivePath(Path.of(SOURCE)) == (TARGET ? FileHelper.asPath(TARGET) : null)

        where:
        SOURCE                      | TARGET
        '/some/file.txt'            | null
        '/data/work/some/file.txt'  | 'http://bucket/data/export/work/some/file.txt'

    }

}
