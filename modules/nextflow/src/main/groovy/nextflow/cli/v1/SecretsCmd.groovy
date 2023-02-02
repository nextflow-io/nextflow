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

package nextflow.cli.v1

import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.cli.CmdSecrets
import nextflow.exception.AbortOperationException

/**
 * CLI `secrets` sub-command (v1)
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
@Parameters(commandDescription = 'Manage pipeline secrets (preview)')
class SecretsCmd extends AbstractCmd implements UsageAware {

    interface SubCmd {
        String getName()
        void apply(List<String> args)
        void usage(List<String> args)
    }

    static public final String NAME = 'secrets'

    private List<SubCmd> commands = (List<SubCmd>)[
        new GetCmd(),
        new SetCmd(),
        new ListCmd(),
        new DeleteCmd()
    ]

    @Parameter(hidden = true)
    List<String> args

    @Override
    String getName() { NAME }

    @Override
    void run() {
        if( !args ) {
            usage()
            return
        }

        final cmd = commands.find { it.name == args[0] }
        if( !cmd ) {
            def matches = commands.collect{ it.name }.closest(args[0])
            def msg = "Unknown secrets sub-command: ${args[0]}"
            if( matches )
                msg += " -- Did you mean one of these?\n" + matches.collect { "  $it"}.join('\n')
            throw new AbortOperationException(msg)
        }

        cmd.apply(args.drop(1))
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

    private void addOption(String fieldName, List<String> args) {
        def annot = this.class.getDeclaredField(fieldName)?.getAnnotation(Parameter)
        if( annot ) {
            args << '  ' + annot.names().join(', ')
            args << '     ' + annot.description()
        }
        else {
            log.debug "Unknown help field: $fieldName"
        }
    }

    @Deprecated
    class PutCmd extends SetCmd {
        @Override
        String getName() { 'put' }

        void apply(List<String> args) {
            log.warn "Put command is deprecated - use 'set' instead'"
            super.apply(args)
        }
    }

    class SetCmd implements SubCmd {
        @Override
        String getName() { 'set' }

        @Override
        void apply(List<String> args) {
            if( args.size() < 1 )
                throw new AbortOperationException("Missing secret name")
            if( args.size() < 2 )
                throw new AbortOperationException("Missing secret value")

            new CmdSecrets().run(CmdSecrets.Command.SET, args)
        }

        @Override
        void usage(List<String> args) {
            args << 'Set a key-pair in the secrets store'
            args << "Usage: nextflow secrets $name <NAME> <VALUE>".toString()
            args << ''
            args << ''
        }
    }

    class GetCmd implements SubCmd {
        @Override
        String getName() { 'get' }

        @Override
        void apply(List<String> args) {
            if( args.size() != 1 )
                throw new AbortOperationException("Wrong number of arguments")

            if( !args[0] )
                throw new AbortOperationException("Missing secret name")

            new CmdSecrets().run(CmdSecrets.Command.GET, args)
        }

        @Override
        void usage(List<String> args) {
            args << 'Get a secret value with the name'
            args << "Usage: nextflow secrets $name <NAME>".toString()
            args << ''
        }
    }

    class ListCmd implements SubCmd {
        @Override
        String getName() { 'list' }

        @Override
        void apply(List<String> args) {
            if( args.size() > 0 )
                throw new AbortOperationException("Wrong number of arguments")

            new CmdSecrets().run(CmdSecrets.Command.LIST, args)
        }

        @Override
        void usage(List<String> args) {
            args << 'List all names in the secrets store'
            args << "Usage: nextflow secrets $name".toString()
            args << ''
        }
    }

    class DeleteCmd implements SubCmd {
        @Override
        String getName() { 'delete' }

        @Override
        void apply(List<String> args) {
            if( args.size() != 1 )
                throw new AbortOperationException("Wrong number of arguments")

            if( !args[0] )
                throw new AbortOperationException("Missing secret name")

            new CmdSecrets().run(CmdSecrets.Command.DELETE, args)
        }

        @Override
        void usage(List<String> args) {
            args << 'Delete an entry from the secrets store'
            args << "Usage: nextflow secrets $name".toString()
            args << ''
            addOption('secretName', args)
            args << ''
        }
    }
}
