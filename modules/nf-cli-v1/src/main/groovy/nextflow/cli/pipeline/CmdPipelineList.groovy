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

import com.beust.jcommander.IParameterValidator
import com.beust.jcommander.Parameter
import com.beust.jcommander.ParameterException
import com.beust.jcommander.Parameters
import groovy.json.JsonOutput
import groovy.transform.CompileStatic
import nextflow.cli.CmdBase
import nextflow.pipeline.PipelineManager
import nextflow.scm.AssetManager
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

    @Parameter(names = ['-g', '-global'], description = 'List the pipelines in the user-level store used by `nextflow run`', arity = 0)
    boolean global

    @TestOnly
    protected Path root

    @Override
    String getName() {
        return 'list'
    }

    @Override
    void run() {
        if( global ) {
            listUserLevel()
            return
        }
        final manager = new PipelineManager(root ?: Path.of('.').toAbsolutePath().normalize(), null)
        final rows = manager.list().collect { p ->
            [
                name   : p.info().name,
                version: p.info().version,
                status : p.info().isModified(p.dir()) ? 'modified' : 'ok'
            ]
        }
        if( output == 'json' ) {
            println JsonOutput.prettyPrint(JsonOutput.toJson(pipelines: rows))
            return
        }
        if( !rows ) {
            println "No pipelines installed"
            final count = AssetManager.list().size()
            if( count )
                println "${count} ${count == 1 ? 'pipeline is' : 'pipelines are'} installed in the user-level store -- use `nextflow pipeline list -g` to list them"
            return
        }
        println "Pipeline".padRight(40) + "Version".padRight(15) + "Status"
        println "-" * 65
        rows.each { r ->
            println "${cell(r.name).padRight(40)}${cell(r.version).padRight(15)}${r.status}"
        }
    }

    private void listUserLevel() {
        final names = AssetManager.list()
        if( output == 'json' ) {
            println JsonOutput.prettyPrint(JsonOutput.toJson(pipelines: names.collect { String it -> [name: it] }))
            return
        }
        if( !names ) {
            println "No pipelines installed in the user-level store"
            return
        }
        names.each { String it -> println it }
    }

    // a damaged .pipeline-info may lack any field
    private static String cell(Object value) {
        return value != null ? value.toString() : '-'
    }
}
