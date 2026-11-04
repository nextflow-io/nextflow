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

package nextflow.cli.pipeline

import java.nio.file.Path
import java.nio.file.Paths

import com.beust.jcommander.IParameterValidator
import com.beust.jcommander.Parameter
import com.beust.jcommander.ParameterException
import com.beust.jcommander.Parameters
import groovy.json.JsonOutput
import groovy.transform.CompileStatic
import nextflow.cli.CmdBase
import nextflow.pipeline.PipelineManager
import nextflow.util.TestOnly

/**
 * Pipeline list subcommand
 */
@CompileStatic
@Parameters(commandDescription = "List installed pipelines")
class CmdPipelineList extends CmdBase {

    @Parameter(names = ['-o', '-output'], description = 'Output mode: text, json', validateWith = OutputModeValidator)
    String output = 'text'

    static class OutputModeValidator implements IParameterValidator {

        private static final List<String> MODES = List.of('text', 'json')

        @Override
        void validate(String name, String value) {
            if( !MODES.contains(value) )
                throw new ParameterException("Output mode must be one of $MODES (found: $value)")
        }
    }

    @TestOnly
    protected Path root

    @Override
    String getName() {
        return 'list'
    }

    @Override
    void run() {
        final manager = new PipelineManager(root ?: Paths.get('.').toAbsolutePath().normalize(), null)
        final rows = manager.list().collect { info ->
            [
                name   : info.name,
                version: info.version,
                source : info.source,
                status : info.isModified(manager.pipelineDir(info.name)) ? 'modified' : 'ok'
            ]
        }
        if( output == 'json' ) {
            println JsonOutput.prettyPrint(JsonOutput.toJson(rows))
            return
        }
        if( !rows ) {
            println "No pipelines installed"
            return
        }
        rows.each { r -> println "${r.name}  ${r.version}  ${r.source}  ${r.status}" }
    }
}
