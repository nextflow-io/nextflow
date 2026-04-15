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

package nextflow.processor

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import nextflow.exception.IllegalArityException
import nextflow.exception.MissingFileException
import nextflow.exception.MissingValueException
import nextflow.script.ScriptType
import nextflow.script.params.CmdEvalParam
import nextflow.script.params.DefaultOutParam
import nextflow.script.params.EnvOutParam
import nextflow.script.params.FileOutParam
import nextflow.script.params.OutParam
import nextflow.script.params.StdOutParam
import nextflow.script.params.TupleOutParam
import nextflow.script.params.ValueOutParam
/**
 * Implements the resolution of task outputs
 * for legacy processes.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class TaskOutputResolverV1 {

    private TaskRun task

    TaskOutputResolverV1(TaskRun task) {
        this.task = task
    }

    void resolve(OutParam param) {
        switch( param ) {
            case StdOutParam:
                collectStdOut((StdOutParam) param)
                break

            case FileOutParam:
                collectOutFiles((FileOutParam) param)
                break

            case ValueOutParam:
                collectOutValue((ValueOutParam) param, task.context)
                break

            case EnvOutParam:
                collectOutEnv((EnvOutParam) param)
                break

            case CmdEvalParam:
                collectOutEval((CmdEvalParam) param)
                break

            case DefaultOutParam:
                task.setOutput(param, DefaultOutParam.Completion.DONE)
                break

            default:
                throw new IllegalArgumentException("Invalid process output: ${param.class.simpleName}")
        }
    }

    /**
     * Resolve a process `env` output.
     *
     * @param param
     */
    protected void collectOutEnv(EnvOutParam param) {
        final value = collectOutEnvMap(task.workDir, null).get(param.name)
        if( value == null && !param.optional )
            throw new MissingValueException("Missing environment variable: ${param.name}")

        task.setOutput(param, value)
    }

    /**
     * Resolve a process `eval` output.
     *
     * @param param
     */
    protected void collectOutEval(CmdEvalParam param) {
        final evalCmds = task.getOutputEvals()
        final value = collectOutEnvMap(task.workDir, evalCmds).get(param.name)
        if( value == null && !param.optional )
            throw new MissingValueException("Missing environment variable: ${param.name}")

        task.setOutput(param, value)
    }

    /**
     * Parse the `.command.env` file which holds the value for `env` and `eval`
     * outputs.
     *
     * @param workDir
     *      The task work directory that contains the `.command.env` file
     * @param evalCmds
     *      A {@link Map} instance containing key-value pairs
     */
    @Memoized(maxCacheSize = 10_000)
    protected Map collectOutEnvMap(Path workDir, Map<String,String> evalCmds) {
        return new TaskEnvCollector(workDir, evalCmds).collect()
    }

    /**
     * Resolve a process `stdout` output.
     *
     * @param param
     */
    protected void collectStdOut(StdOutParam param) {
        final stdout = task.@stdout

        if( stdout == null && task.type == ScriptType.SCRIPTLET )
            throw new IllegalArgumentException("Missing 'stdout' for process > ${task.lazyName()}")

        if( stdout instanceof Path && !stdout.exists() )
            throw new MissingFileException("Missing 'stdout' file: ${stdout.toUriString()} for process > ${task.lazyName()}")

        final result = stdout instanceof Path ? stdout.text : stdout?.toString()
        task.setOutput(param, result)
    }

    /**
     * Resolve a process `file` or `path` output.
     *
     * @param param
     */
    protected void collectOutFiles(FileOutParam param) {

        // `file` outputs can specify multiple file patterns separated by `:`
        final filePatterns = param.getFilePatterns(task.context, task.workDir)
        final opts = [
            followLinks: param.followLinks,
            glob: param.glob,
            hidden: param.hidden,
            includeInputs: param.includeInputs,
            maxDepth: param.maxDepth,
            optional: param.optional || param.arity?.min == 0,
            type: param.type,
        ]
        final allFiles = collectOutFiles0(filePatterns, opts)

        if( !param.isValidArity(allFiles.size()) )
            throw new IllegalArityException("Incorrect number of output files for process `${task.lazyName()}` -- expected ${param.arity}, found ${allFiles.size()}")

        final result = allFiles.size() == 1 && param.isSingle() ? allFiles[0] : allFiles
        task.setOutput(param, result)
    }

    protected List<Path> collectOutFiles0(List<String> filePatterns, Map<String,?> opts) {
        return new TaskFileCollector(filePatterns, opts, task).collect()
    }

    /**
     * Resolve a process `val` output.
     *
     * @param param
     * @param ctx
     */
    protected void collectOutValue(ValueOutParam param, Map ctx) {
        try {
            task.setOutput(param, param.resolve(ctx))
        }
        catch( MissingPropertyException e ) {
            throw new MissingValueException("Missing value declared as output parameter: ${e.property}")
        }
    }

}
