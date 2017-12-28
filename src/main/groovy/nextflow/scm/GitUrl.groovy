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

import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString

/**
 * Models Git URL
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@ToString(includeFields=true, includes='protocol,user,domain,project')
@EqualsAndHashCode(includeFields=true, includes='protocol,user,domain,project')
class GitUrl {

    private String protocol

    private String user

    private String domain

    private String project

    /**
     * Creates a git url object specifying the url components
     *
     * @param protocol The url protocol i.e. {@code http}
     * @param user The remote user name
     * @param domain Either the server domain {@code http} or the path location for local file URL {@code /usr/local/repo}
     * @param project The git project name i.e. {@code foo/bar}
     */
    GitUrl( String protocol, String user, String domain, String project ) {
        this.protocol = protocol
        this.user = user
        this.domain = domain
        this.project = project
    }

    /**
     * Creates a Git URL object parsing the URL string
     *
     * @param url The git URL to parse i.e. {@code https://github.com/cbcrg/grape.git}
     */
    GitUrl( String url ) {
        if( !url ) throw new IllegalArgumentException("Git URL cannot be empty")

        if( !url.endsWith('.git') )
            url += '.git'

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
            this.domain = url.substring(p+1)
            user = url.substring(0,p)
        }
        else {
            this.domain = url
        }
    }

    private parse1( String url, int p ) {
        url = url.substring(p+1)
        p = url.indexOf(':')
        if( p != -1 ) {
            project = url.substring(p+1)
            if( project.endsWith('.git'))
                project = project.substring(0,project.size()-4)
            this.domain = url.substring(0,p)
        }
        else {
            this.domain = url
        }

    }


    private parseFile( String str ) {
        if( str.endsWith('.git') ) {
            str = str.substring(0, str.size()-4)
            def items = str.tokenize('/')
            project = items[-1]
            if( items.size()>1 ) {
                this.domain = '/' + items[0..-2].join('/')
            }
        }
        else {
            this.domain = str
        }
    }

    /**
     * @return
     *      The url protocol e.g. http, https, git, ssh or file. Not it fallbacks on
     *      @{code file} when an absolute path is specified without the {@coe file:} prefix
     */
    String getProtocol() { protocol }

    /**
     * @return The user name eventually specified in the url string
     */
    String getUser() { user }

    /**
     * @return
     *      The URL domain e.g. given the URL {@code https://gitlab.com/pditommaso/hello.git} the string {@code gitlab.com}
     *      is return. When the url is file system path the path withoit the last component is returned
     *      e.g. given {@code file:///opt/git/project.git} the following string is returned {@code /opt/git/}
     */
    String getDomain() { this.domain }

    /**
     * @return
     *      The repository qualified name e.g. given the project repository {@code https://github.com/cbcrg/piper-nf.git}
     *      returns the string {@code cbcrg/piper-nf}.
     */
    String getProject() { project }

}
