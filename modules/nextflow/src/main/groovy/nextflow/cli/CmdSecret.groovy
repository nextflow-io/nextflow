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
 *
 */

package nextflow.cli

import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.exception.AbortOperationException
import nextflow.plugin.Plugins
import nextflow.secret.SecretsLoader
import nextflow.secret.SecretsProvider
/**
 * Implements the {@code secret} command
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class CmdSecret {

    static public final String NAME = 'secrets'

    enum Command {
        GET,
        SET,
        LIST,
        DELETE
    }

    @Parameters(commandDescription = "Manage pipeline secrets")
    static class V1 extends CmdBase implements UsageAware {

        interface SubCmd {
            String getName()
            void apply(List<String> result)
            void usage(List<String> result)
        }

        private List<SubCmd> commands = []

        String getName() {
            return NAME
        }

        @Parameter(hidden = true)
        List<String> args = []

        V1() {
            commands.add( new GetCmd() )
            commands.add( new PutCmd() )
            commands.add( new SetCmd() )
            commands.add( new ListCmd() )
            commands.add( new DeleteCmd() )
        }

        /**
         * Print the command usage help
         */
        @Override
        void usage() {
            usage(args)
        }

        /**
         * Print the command usage help
         *
         * @param args The arguments as entered by the user
         */
        @Override
        void usage(List<String> args) {

            List<String> result = []
            if( !args ) {
                result << this.getClass().getAnnotation(Parameters).commandDescription()
                result << 'Usage: nextflow secrets <sub-command> [options]'
                result << ''
                result << 'Commands:'
                commands.collect{ it.name }.sort().each { result << "  $it".toString()  }
                result << ''
            }
            else {
                def sub = commands.find { it.name == args[0] }
                if( sub )
                    sub.usage(result)
                else {
                    throw new AbortOperationException("Unknown secrets sub-command: ${args[0]}")
                }
            }

            println result.join('\n').toString()
        }

        @Override
        void run() {
            if( !args ) {
                usage()
                return
            }

            getCmd(args).apply(args.drop(1))
        }

        protected SubCmd getCmd(List<String> args) {

            def cmd = commands.find { it.name == args[0] }
            if( cmd ) {
                return cmd
            }

            def matches = commands.collect{ it.name }.closest(args[0])
            def msg = "Unknown secrets sub-command: ${args[0]}"
            if( matches )
                msg += " -- Did you mean one of these?\n" + matches.collect { "  $it"}.join('\n')
            throw new AbortOperationException(msg)
        }

        private void addOption(String fieldName, List<String> result) {
            def annot = this.class.getDeclaredField(fieldName)?.getAnnotation(Parameter)
            if( annot ) {
                result << '  ' + annot.names().join(', ')
                result << '     ' + annot.description()
            }
            else {
                log.debug "Unknown help field: $fieldName"
            }
        }

        @Deprecated
        class PutCmd extends SetCmd {
            @Override
            String getName() { 'put' }

            @Override
            void apply(List<String> result) {
                log.warn "Put command is deprecated - use 'set' instead'"
                super.apply(result)
            }
        }

        class SetCmd implements SubCmd {

            @Override
            String getName() { 'set' }

            @Override
            void apply(List<String> result) {
                if( result.size() < 1 )
                    throw new AbortOperationException("Missing secret name")
                if( result.size() < 2 )
                    throw new AbortOperationException("Missing secret value")

                new CmdSecret().run(Command.SET, result)
            }

            @Override
            void usage(List<String> result) {
                result << 'Set a key-pair in the secrets store'
                result << "Usage: nextflow secrets $name <NAME> <VALUE>".toString()
                result << ''
                result << ''
            }
        }

        class GetCmd implements SubCmd {

            @Override
            String getName() { 'get' }

            @Override
            void apply(List<String> result) {
                if( result.size() != 1 )
                    throw new AbortOperationException("Wrong number of arguments")

                if( !result.first() )
                    throw new AbortOperationException("Missing secret name")

                new CmdSecret().run(Command.GET, result)
            }

            @Override
            void usage(List<String> result) {
                result << 'Get a secret value with the name'
                result << "Usage: nextflow secrets $name <NAME>".toString()
                result << ''
            }
        }

        class ListCmd implements SubCmd {
            @Override
            String getName() { 'list' }

            @Override
            void apply(List<String> result) {
                if( result.size() > 0 )
                    throw new AbortOperationException("Wrong number of arguments")

                new CmdSecret().run(Command.LIST, result)
            }

            @Override
            void usage(List<String> result) {
                result << 'List all names in the secrets store'
                result << "Usage: nextflow secrets $name".toString()
                result << ''
            }
        }

        class DeleteCmd implements SubCmd {
            @Override
            String getName() { 'delete' }

            @Override
            void apply(List<String> result) {
                if( result.size() != 1 )
                    throw new AbortOperationException("Wrong number of arguments")

                if( !result.first() )
                    throw new AbortOperationException("Missing secret name")

                new CmdSecret().run(Command.DELETE, result)
            }

            @Override
            void usage(List<String> result) {
                result << 'Delete an entry from the secrets store'
                result << "Usage: nextflow secrets $name".toString()
                result << ''
                addOption('secretName', result)
                result << ''
            }
        }
    }

    private SecretsProvider provider

    void run(Command command, List<String> args) {
        // setup the plugins system and load the secrets provider
        Plugins.init()
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
