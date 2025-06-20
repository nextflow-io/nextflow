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

package nextflow.scm

import groovy.transform.CompileStatic
/**
 * Options for interacting with a git repository provider i.e. GitHub or BitBucket
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
record HubOptions(
    String provider,
    String user
) {

    String getPassword() {
        if( !user )
            return null

        final p = user.indexOf(':')
        if( p != -1 )
            return user.substring(p + 1)

        final console = System.console()
        if( !console )
            return null

        print "Enter your $provider password: "
        return new String(console.readPassword())
    }

    String getUser() {
        if( !user )
            return user

        final p = user.indexOf(':')
        return p != -1 ? user.substring(0, p) : user
    }
}
