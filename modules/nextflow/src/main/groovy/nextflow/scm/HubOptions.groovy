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

package nextflow.scm

import groovy.transform.CompileStatic
/**
 * Options for interacting with a git repository provider i.e. GitHub or BitBucket
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class HubOptions {

    private final String provider

    private final String credentials

    HubOptions(String provider, String credentials) {
        this.provider = provider
        this.credentials = credentials
    }

    String getProvider() { provider }

    String getUser() {
        if( !credentials )
            return credentials

        final p = credentials.indexOf(':')
        return p != -1 ? credentials.substring(0, p) : credentials
    }

    String getPassword() {
        if( !credentials )
            return null

        final p = credentials.indexOf(':')
        if( p != -1 )
            return credentials.substring(p + 1)

        final console = System.console()
        if( !console )
            return null

        print "Enter your $provider password: "
        return new String(console.readPassword())
    }

    /**
     * @return A copy of these options bound to the given provider, used to resolve the
     * actual hub provider name (e.g. for the password prompt) while preserving the raw
     * credentials. Keeping {@code credentials} private avoids exposing the raw
     * {@code user:password} string as a public accessor.
     */
    HubOptions withProvider(String provider) {
        new HubOptions(provider, credentials)
    }
}
