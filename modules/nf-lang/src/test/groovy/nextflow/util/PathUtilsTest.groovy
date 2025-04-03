/*
 * Copyright 2024-2025, Seqera Labs
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

import java.nio.file.Path

import spock.lang.Specification

/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class PathUtilsTest extends Specification {

    def 'should determine whether a path is excluded' () {
        given:
        def patterns = [ '.git', '.nf-test', 'work', 'foo/bar' ]

        expect:
        result == PathUtils.isExcluded(Path.of(path), patterns)

        where:
        path                            | result
        '.git/main.nf'                  | true
        '.nf-test/mock.nf'              | true
        'modules/foo/.nf-test/mock.nf'  | true
        'modules/foo/bar/main.nf'       | true
        'work/01/234567/main.nf'        | true
        'main.nf'                       | false
        'modules/foo/main.nf'           | false
        'nextflow.config'               | false
        'conf/base.config'              | false
    }

}
