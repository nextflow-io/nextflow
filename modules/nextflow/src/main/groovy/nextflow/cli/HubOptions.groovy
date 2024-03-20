/*
 * Copyright 2013-2024, Seqera Labs
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

package nextflow.cli

import com.beust.jcommander.Parameter
import groovy.transform.CompileStatic
/**
  * Defines the command line parameters for command that need to interact with a pipeline service hub i.e. GitHub or BitBucket
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@CompileStatic
trait HubOptions {

    abstract String getHubProvider()

    abstract String getHubUserCli()

    abstract void setHubProvider(String hub)

    static trait V1 implements HubOptions {

        @Parameter(names=['-hub'], description = "Service hub where the project is hosted")
        String hubProvider

        @Parameter(names='-user', description = 'Private repository user name')
        String hubUserCli

    }

    /**
     * Return the password provided on the command line or stop allowing the user to enter it on the console
     *
     * @return The password entered or {@code null} if no user has been entered
     */
    String getHubPassword() {

        if( !hubUserCli )
            return null

        def p = hubUserCli.indexOf(':')
        if( p != -1 )
            return hubUserCli.substring(p+1)

        def console = System.console()
        if( !console )
            return null

        print "Enter your $hubProvider password: "
        char[] pwd = console.readPassword()
        new String(pwd)
    }

    String getHubUser() {
        if( !hubUserCli ) {
            return hubUserCli
        }

        def p = hubUserCli.indexOf(':')
        return p != -1 ? hubUserCli.substring(0,p) : hubUserCli
    }

}
