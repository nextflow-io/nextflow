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
import nextflow.cli.pipeline.CmdPipelineInstall
import nextflow.cli.pipeline.CmdPipelineList
import nextflow.cli.pipeline.CmdPipelineRemove
import nextflow.cli.pipeline.CmdPipelineSearch
import nextflow.cli.pipeline.CmdPipelineView

/**
 * Implements `pipeline` command
 */
@CompileStatic
@Parameters(commandDescription = "Manage registry pipelines")
class CmdPipeline extends CmdGroup {

    static final public String NAME = 'pipeline'

    static final List<CmdBase> commands = new ArrayList<>()

    static {
        commands << new CmdPipelineInstall()
        commands << new CmdPipelineList()
        commands << new CmdPipelineRemove()
        commands << new CmdPipelineSearch()
        commands << new CmdPipelineView()
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
