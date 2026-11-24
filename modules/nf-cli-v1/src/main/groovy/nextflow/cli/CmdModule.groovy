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

package nextflow.cli

import com.beust.jcommander.Parameters
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.cli.module.CmdModuleCreate
import nextflow.cli.module.CmdModuleView
import nextflow.cli.module.CmdModuleInstall
import nextflow.cli.module.CmdModuleList
import nextflow.cli.module.CmdModulePublish
import nextflow.cli.module.CmdModuleRemove
import nextflow.cli.module.CmdModuleRun
import nextflow.cli.module.CmdModuleSearch
import nextflow.cli.module.CmdModuleSpec
import nextflow.cli.module.CmdModuleValidate

/**
 * Implements `module` command
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@CompileStatic
@Slf4j
@Parameters(commandDescription = "Manage Nextflow modules", commandNames = ['modules'])
class CmdModule extends CmdGroup {

    static final public String NAME = 'module'

    static final List<CmdBase> commands = new ArrayList<>()

    static {
        commands << new CmdModuleCreate()
        commands << new CmdModuleInstall()
        commands << new CmdModuleRun()
        commands << new CmdModuleList()
        commands << new CmdModuleRemove()
        commands << new CmdModuleSearch()
        commands << new CmdModuleView()
        commands << new CmdModulePublish()
        commands << new CmdModuleSpec()
        commands << new CmdModuleValidate()
    }

    @Override
    protected List<CmdBase> subcommands() {
        return commands
    }

    @Override
    String getName() {
        return NAME
    }
}
