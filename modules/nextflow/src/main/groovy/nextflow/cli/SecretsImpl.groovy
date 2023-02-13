/*
 * Copyright 2020-2022, Seqera Labs
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
 *
 */

package nextflow.cli

import groovy.transform.CompileStatic
import nextflow.plugin.Plugins
import nextflow.secret.SecretsLoader
import nextflow.secret.SecretsProvider

/**
 * CLI `secrets` sub-command
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class SecretsImpl {

    enum Command {
        GET,
        SET,
        LIST,
        DELETE
    }

    private SecretsProvider provider

    void run(Command command, List<String> args) {
        // setup the plugins system and load the secrets provider
        Plugins.setup()
        provider = SecretsLoader.instance.load()

        // run the command
        try {
            switch( command ) {
                case GET:
                    get(args[0])
                    break
                case SET:
                    set(args[0], args[1])
                    break
                case LIST:
                    list()
                    break
                case DELETE:
                    delete(args[0])
                    break
            }
        }
        finally {
            // close the provider
            provider?.close()
        }
    }

    void get(String name) {
        println provider.getSecret(name)?.value
    }

    void set(String name, String value) {
        provider.putSecret(name, value)
    }

    void list() {
        final names = new ArrayList(provider.listSecretsNames()).sort()
        if( names.size() == 0 ) {
            println "no secrets available"
        }

        for( String it : names ) {
            println it
        }
    }

    void delete(String name) {
        provider.removeSecret(name)
    }
}
