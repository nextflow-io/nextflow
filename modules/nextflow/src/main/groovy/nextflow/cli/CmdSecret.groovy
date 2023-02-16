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
import groovy.util.logging.Slf4j
import nextflow.exception.AbortOperationException
import nextflow.plugin.Plugins
import nextflow.secret.SecretsLoader
import nextflow.secret.SecretsProvider
import picocli.CommandLine.Command
import picocli.CommandLine.Parameters
/**
 * Implements the {@code secrets} command
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
@Command(name = 'secrets', description = "Manage pipeline secrets (preview)")
class CmdSecret extends CmdBase {

    enum SubCommand {
        GET,
        SET,
        LIST,
        DELETE
    }

    @Command(description = 'Set a key-pair in the secrets store')
    void get(
            @Parameters(paramLabel = '<name>') String name) {
        run0(SubCommand.GET, [ name ])
    }

    @Command(description = 'Get a secret value with the name')
    void set(
            @Parameters(paramLabel = '<name>') String name,
            @Parameters(paramLabel = '<value>') String value) {
        run0(SubCommand.SET, [ name, value ])
    }

    @Command(description = 'List all names in the secrets store')
    void list() {
        run0(SubCommand.LIST, [])
    }

    @Command(description = 'Delete an entry from the secrets store')
    void delete(
            @Parameters(paramLabel = '<name>') String name) {
        run0(SubCommand.DELETE, [ name ])
    }

    private SecretsProvider provider

    private void run0(SubCommand command, List<String> args) {
        // setup the plugins system and load the secrets provider
        Plugins.setup()
        provider = SecretsLoader.instance.load()

        // run the command
        try {
            switch( command ) {
                case GET:
                    get0(args[0])
                    break
                case SET:
                    set0(args[0], args[1])
                    break
                case LIST:
                    list0()
                    break
                case DELETE:
                    delete0(args[0])
                    break
            }
        }
        finally {
            // close the provider
            provider?.close()
        }
    }

    private void get0(String name) {
        println provider.getSecret(name)?.value
    }

    private void set0(String name, String value) {
        provider.putSecret(name, value)
    }

    private void list0() {
        final names = new ArrayList(provider.listSecretsNames()).sort()
        if( names.size() == 0 ) {
            println "no secrets available"
        }

        for( String it : names ) {
            println it
        }
    }

    private void delete0(String name) {
        provider.removeSecret(name)
    }
}
