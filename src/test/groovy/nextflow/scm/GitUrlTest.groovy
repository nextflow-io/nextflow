/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
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
        result.project == project
        result.domain == domain
        result.user == user
        result == new GitUrl(protocol, user, domain, project)

        where:
        url                                             | protocol  | user  | domain            | project
        'https://example.com/gitproject.git'            | 'https'   | null  | 'example.com'     | 'gitproject'
        'https://gitlab.com/pditommaso/hello.git'       | 'https'   | null  | 'gitlab.com'      | 'pditommaso/hello'
        'http://github.com/pditommaso/hello.git'        | 'http'    | null  | 'github.com'      | 'pditommaso/hello'
        'https://yo@github.com/pditommaso/hola.git'     | 'https'   | 'yo'  | 'github.com'      | 'pditommaso/hola'
        'git@gitlab.com:pditommaso/hello.git'           | 'git'     | null  | 'gitlab.com'      | 'pditommaso/hello'
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
