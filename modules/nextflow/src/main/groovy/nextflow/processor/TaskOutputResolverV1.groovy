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
import nextflow.script.params.BaseOutParam
import nextflow.script.params.CmdEvalParam
import nextflow.script.params.DefaultOutParam
import nextflow.script.params.EnvOutParam
import nextflow.script.params.FileOutParam
import nextflow.script.params.OptionalParam
import nextflow.script.params.OutParam
import nextflow.script.params.StdOutParam
import nextflow.script.params.ValueOutParam
/**
 * Implements the resolution of task outputs for legacy processes.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class TaskOutputResolverV1 {

    private TaskRun task

    /**
     * The directory in which the task outputs are collected, i.e.
     * the store directory when specified, otherwise the task
     * work directory.
     */
    private Path workDir

    TaskOutputResolverV1(TaskRun task) {
        this.task = task
        this.workDir = task.getTargetDir()
    }

    /**
     * Resolve an output param and set the resulting value
     * into the task.
     *
     * @param param
     */
    void resolve(OutParam param) {

        switch( param ) {
            case StdOutParam:
                collectStdOut((StdOutParam)param)
                break

            case FileOutParam:
                collectOutFiles((FileOutParam)param)
                break

            case ValueOutParam:
                collectOutValue((ValueOutParam)param)
                break

            case EnvOutParam:
            case CmdEvalParam:
                collectOutEnv((BaseOutParam)param)
                break

            case DefaultOutParam:
                task.setOutput(param, DefaultOutParam.Completion.DONE)
                break

            default:
                throw new IllegalArgumentException("Illegal output parameter: ${param.class.simpleName}")
        }
    }

    /**
     * Collects the process 'std output'
     *
     * @param param
     */
    protected void collectStdOut(StdOutParam param) {
        final stdout = task.@stdout

        if( stdout == null && task.type == ScriptType.SCRIPTLET )
            throw new IllegalArgumentException("Missing 'stdout' for process > ${task.lazyName()}")

        if( stdout instanceof Path && !stdout.exists() )
            throw new MissingFileException("Missing 'stdout' file: ${stdout.toUriString()} for process > ${task.lazyName()}")

        task.setOutput(param, stdout)
    }

    /**
     * Collects a process `file` or `path` output.
     *
     * @param param
     */
    protected void collectOutFiles(FileOutParam param) {

        // type file parameter can contain a multiple files pattern separating them with a special character
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

        task.setOutput( param, allFiles.size()==1 && param.isSingle() ? allFiles[0] : allFiles )
    }

    protected List<Path> collectOutFiles0(List<String> filePatterns, Map<String,?> opts) {
        return new TaskFileCollector(filePatterns, opts, task).collect()
    }

    /**
     * Collects a process `val` output.
     *
     * @param param
     */
    protected void collectOutValue(ValueOutParam param) {

        try {
            // fetch the output value
            final val = param.resolve(task.context)
            // set into the output set
            task.setOutput(param, val)
            // trace the result
            log.trace "Collecting param: ${param.name}; value: ${val}"
        }
        catch( MissingPropertyException e ) {
            throw new MissingValueException("Missing value declared as output parameter: ${e.property}")
        }
    }

    /**
     * Collects a process `env` or `eval` output.
     *
     * @param param
     */
    protected void collectOutEnv(BaseOutParam param) {

        // fetch the output value
        final outCmds = param instanceof CmdEvalParam ? task.getOutputEvals() : null
        final val = collectOutEnvMap(outCmds).get(param.name)
        if( val == null && !((OptionalParam)param).optional )
            throw new MissingValueException("Missing environment variable: $param.name")
        // set into the output set
        task.setOutput(param, val)
        // trace the result
        log.trace "Collecting param: ${param.name}; value: ${val}"
    }

    /**
     * Parse the `.command.env` file which holds the value for `env` and `eval`
     * output types
     *
     * @param outEvals
     *      A {@link Map} instance containing key-value pairs
     */
    @Memoized
    protected Map collectOutEnvMap(Map<String,String> outEvals) {
        return new TaskEnvCollector(workDir, outEvals).collect()
    }

}
