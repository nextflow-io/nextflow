/*
 * Copyright (c) 2013-2014, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2014, Paolo Di Tommaso and the respective authors.
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

import com.beust.jcommander.IValueValidator
import com.beust.jcommander.Parameter
import com.beust.jcommander.ParameterException
import groovy.transform.CompileStatic
import nextflow.scm.AssetManager
/**
  * Defines the command line parameters for command that need to interact with a pipeline service hub i.e. GitHub or BitBucket
  *
  * Author Maria Chatzou
  * Author Paolo Di Tommaso
  */
@CompileStatic
trait HubParams {

    @Parameter(names=['-hub'], description = "Pipeline service hub where the code is hosted - It can be either 'github' or 'bitbucket'", validateValueWith = HubValidator)
    String hub_provider = AssetManager.DEFAULT_HUB

    @Parameter(names='-user', description = 'Private repository user name')
    String hub_user

    @Parameter(names='-pwd', description = 'Private repository password')
    String hub_pwd

    /**
     * Return the password provided on the command line or stop allowing the user to enter it on the console
     *
     * @return The password entered or {@code null} if no user has been entered
     */
    String getHubPassword() {

        if( !hub_user ) {
            return null
        }

        if( hub_pwd ) {
            return hub_pwd
        }

        print "Enter your $hub_provider password: "
        char[] pwd = System.console().readPassword()
        new String(pwd)
    }


    @CompileStatic
    static class HubValidator implements IValueValidator<String> {

        final static List VALID_HUB = ['github','bitbucket']

        @Override
        void validate(String name, String value) throws ParameterException {
            if( name != '-hub' ) return
            if( value in VALID_HUB ) return
            throw new ParameterException("Not a valid '-hub' provider -- It must be either: ${VALID_HUB.join(', ')}")
        }

    }

}