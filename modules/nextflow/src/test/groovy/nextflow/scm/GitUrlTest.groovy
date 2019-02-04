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

package nextflow.scm

import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class GitUrlTest extends Specification {

    @Unroll
    def 'should parse github url: #url' () {

        when:
        def result = new GitUrl(url)
        then:
        result.protocol == protocol
        result.path == path
        result.domain == domain
        result.user == user
        result == new GitUrl(protocol, user, domain, path)

        where:
        url                                             | protocol  | user  | domain            | path
        'https://example.com/gitproject.git'            | 'https'   | null  | 'example.com'     | 'gitproject'
        'https://gitlab.com/pditommaso/hello.git'       | 'https'   | null  | 'gitlab.com'      | 'pditommaso/hello'
        'http://github.com/pditommaso/hello.git'        | 'http'    | null  | 'github.com'      | 'pditommaso/hello'
        'https://gitlab.com/sub1/pditommaso/hello.git'  | 'https'   | null  | 'gitlab.com'      | 'sub1/pditommaso/hello'
        'https://yo@github.com/pditommaso/hola.git'     | 'https'   | 'yo'  | 'github.com'      | 'pditommaso/hola'
        'git@gitlab.com:pditommaso/hello.git'           | 'git'     | null  | 'gitlab.com'      | 'pditommaso/hello'
        'git@gitlab.com:dir1/pditommaso/hello.git'      | 'git'     | null  | 'gitlab.com'      | 'dir1/pditommaso/hello'
        'ssh://me@server/project.git'                   | 'ssh'     | 'me'  | 'server'          | 'project'
        'me@server:project.git'                         | 'ssh'     | 'me'  | 'server'          | 'project'
        'file:///opt/git/project.git'                   | 'file'    | null  | '/opt/git'        | 'project'
        'file:/opt/git/project.git'                     | 'file'    | null  | '/opt/git'        | 'project'
        '/opt/git/project.git'                          | 'file'    | null  | '/opt/git'        | 'project'
        '/opt/git/project/.git'                         | 'file'    | null  | '/opt/git'        | 'project'

    }

    def 'should throw IllegalArgumentException' () {
        when:
        new GitUrl('foo.com/bar')
        then:
        thrown(IllegalArgumentException)
    }


}
