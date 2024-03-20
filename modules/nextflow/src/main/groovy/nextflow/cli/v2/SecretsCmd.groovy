/*
 * Copyright 2013-2023, Seqera Labs
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

package nextflow.cli.v2

import groovy.transform.CompileStatic
import nextflow.cli.CmdSecret
import picocli.CommandLine.Command
import picocli.CommandLine.Parameters

/**
 * CLI `secrets` sub-command (v2)
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
@Command(
    name = 'secrets',
    description = 'Manage pipeline secrets (preview)'
)
class SecretsCmd extends AbstractCmd {

    @Command(description = 'Set a key-pair in the secrets store')
    void get(
            @Parameters(paramLabel = '<name>') String name) {
        new CmdSecret().run(CmdSecret.Command.GET, [ name ])
    }

    @Command(description = 'Get a secret value with the name')
    void set(
            @Parameters(paramLabel = '<name>') String name,
            @Parameters(paramLabel = '<value>') String value) {
        new CmdSecret().run(CmdSecret.Command.SET, [ name, value ])
    }

    @Command(description = 'List all names in the secrets store')
    void list() {
        new CmdSecret().run(CmdSecret.Command.LIST, [])
    }

    @Command(description = 'Delete an entry from the secrets store')
    void delete(
            @Parameters(paramLabel = '<name>') String name) {
        new CmdSecret().run(CmdSecret.Command.DELETE, [ name ])
    }

}
