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

package nextflow.util

import java.nio.file.Files
import java.nio.file.Path

import spock.lang.Specification

import static test.TestUtils.deleteDir
import static test.TestUtils.tempDirInMem
import static test.TestUtils.tempFile

/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class PathUtilsTest extends Specification {

    def 'should visit all files in a directory that satisfy a condition' () {
        given:
        def root = tempDirInMem()
        def files = [
            'main.nf',
            'nextflow.config',
            'modules/index/main.nf',
            'modules/fastqc/main.nf',
            'modules/quant/main.nf',
            'modules/multiqc/main.nf',
            'modules/rnaseq.nf'
        ]
        files.each { filename ->
            tempFile(root, filename)
        }

        when:
        def result = [] as Set
        PathUtils.visitFiles(
            root,
            path -> Files.isDirectory(path) || path.toString().endsWith('.nf'),
            path -> result.add(root.relativize(path).toString()) )
        then:
        result == files.findAll { filename -> filename.endsWith('.nf') }.toSet()

        cleanup:
        deleteDir(root)
    }

    def 'should determine whether a path is excluded' () {
        given:
        def patterns = [ '.git', '.nf-test', 'work', 'foo/bar', 'nf-test.config' ]

        expect:
        result == PathUtils.isExcluded(Path.of(path), patterns)

        where:
        path                            | result
        '.git'                          | true
        '.nf-test'                      | true
        'modules/foo/.nf-test'          | true
        'modules/foo/bar'               | true
        'work'                          | true
        'nf-test.config'                | true
        'modules/foo/nf-test.config'    | true
        'main.nf'                       | false
        'modules/foo/main.nf'           | false
        'subworkflows/foobar/main.nf'   | false
        'nextflow.config'               | false
        'conf/base.config'              | false
    }

}
