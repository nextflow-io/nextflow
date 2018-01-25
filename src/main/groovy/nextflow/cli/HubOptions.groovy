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

package nextflow.cli

import com.beust.jcommander.Parameter
import groovy.transform.CompileStatic
import picocli.CommandLine

/**
  * Defines the command line parameters for command that need to interact with a pipeline service hub i.e. GitHub or BitBucket
  *
  * @author Maria Chatzou
  * @author Paolo Di Tommaso
  */

@CompileStatic
//TODO needs @CommandLine.Command ?
trait HubOptions {

    //@Parameter(names=['-hub'], description = "Service hub where the project is hosted")
    @CommandLine.Option(names = ['--hub'], description = "Service hub where the project is hosted")
    String hubProvider

    //@Parameter(names='-user', description = 'Private repository user name')
    @CommandLine.Option(names=['--user'], description = 'Private repository user name')
    String hubUser

    /**
     * Return the password provided on the command line or stop allowing the user to enter it on the console
     *
     * @return The password entered or {@code null} if no user has been entered
     */
    String getHubPassword() {

        if( !hubUser )
            return null

        def p = hubUser.indexOf(':')
        if( p != -1 )
            return hubUser.substring(p+1)

        def console = System.console()
        if( !console )
            return null

        print "Enter your $hubProvider password: "
        char[] pwd = console.readPassword()
        new String(pwd)
    }

    String getHubUser() {
        if(!hubUser) {
            return hubUser
        }

        def p = hubUser.indexOf(':')
        return p != -1 ? hubUser.substring(0,p) : hubUser
    }

}