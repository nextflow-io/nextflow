/*
 * Copyright (c) 2013-2015, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2015, Paolo Di Tommaso and the respective authors.
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

import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString
/**
 * Parse a git URL
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@ToString()
@EqualsAndHashCode(includeFields=true, includes='protocol,user,location,project')
class GitUrlParser {

    private String protocol

    private String user

    private String location

    private String project

    GitUrlParser( String protocol, String user, String location, String project ) {
        this.protocol = protocol
        this.user = user
        this.location = location
        this.project = project
    }

    GitUrlParser( String url ) {

        def p = url.indexOf('://')
        if( p != -1 ) {
            parse0(url, p)
            return
        }

        else if( ( p = url.indexOf('@')) != -1) {
            user = url.substring(0,p)
            if( user == 'git' ) {
                protocol = 'git'
                user = null
            }
            else {
                protocol = 'ssh'
            }
            parse1(url, p)
            return
        }

        else if( url.startsWith('file:') ) {
            protocol = 'file'
            parseFile(url.substring(5))
            return
        }

        else if( url.startsWith('/') ) {
            protocol = 'file'
            parseFile(url)
            return
        }

        throw new IllegalArgumentException("Malformed repository url: $url")
    }


    private parse0( String url, int p ) {
        // split protocol and url
        protocol = url.substring(0,p)

        url = url.substring(p+3)
        if( protocol == 'file' ) {
            parseFile(url)
            return
        }

        // remove the remaining part
        p = url.indexOf('/')
        if( p != -1 ) {
            project = url.substring(p)
            if( project.endsWith('.git'))
                project = project.substring(0,project.size()-4)
            if( project.startsWith('/') )
                project = project.substring(1)
            url = url.substring(0,p)
        }
        // remove the user name
        p = url.indexOf('@')
        if( p != -1 ) {
            location = url.substring(p+1)
            user = url.substring(0,p)
        }
        else {
            location = url
        }
    }

    private parse1( String url, int p ) {
        url = url.substring(p+1)
        p = url.indexOf(':')
        if( p != -1 ) {
            project = url.substring(p+1)
            if( project.endsWith('.git'))
                project = project.substring(0,project.size()-4)
            location = url.substring(0,p)
        }
        else {
            location = url
        }

    }


    private parseFile( String str ) {
        if( str.endsWith('.git') ) {
            str = str.substring(0, str.size()-4)
            def items = str.tokenize('/')
            project = items[-1]
            if( items.size()>1 ) {
                location = '/' + items[0..-2].join('/')
            }
        }
        else {
            location = str
        }
    }


    String getProtocol() { protocol }

    String getUser() { user }

    String getLocation() { location }

    String getProject() { project }

}
