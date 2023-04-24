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
 */

package nextflow.executor

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.exception.ProcessFailedException
import nextflow.exception.ProcessNonZeroExitStatusException
import nextflow.processor.TaskHandler
import nextflow.processor.TaskStatus
import nextflow.util.CmdLineHelper
/**
 * Handles the execution of an array job for any grid executor.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class ArrayGridTaskHandler extends ArrayTaskHandler implements SubmitRetryAware {

    final AbstractGridExecutor executor

    private jobId

    ArrayGridTaskHandler(List<TaskHandler> array, AbstractGridExecutor executor) {
        super(array)

        this.executor = executor
    }

    @Override
    void submit() {
        ProcessBuilder builder = null
        try {
            // -- create the array job script
            final arrayScript = executor.createArrayTaskWrapper(this)

            // -- create the submit command
            builder = createProcessBuilder()

            // -- submit the array job with a retryable strategy
            final result = safeExecute( () -> processStart(builder, arrayScript) )

            // -- save the job id
            this.setJobId(executor.parseJobId(result))
            this.setStatus(TaskStatus.SUBMITTED)

            log.debug "[${executor.name.toUpperCase()}] submitted array job > jobId: ${jobId}"
        }
        catch( Exception e ) {
            // update task exit status and message
            if( e instanceof ProcessNonZeroExitStatusException ) {
                for( TaskHandler handler : array ) {
                    handler.task.exitStatus = e.getExitStatus()
                    handler.task.stdout = e.getReason()
                    handler.task.script = e.getCommand()
                }
            }
            else {
                for( TaskHandler handler : array )
                    handler.task.script = builder ? CmdLineHelper.toLine(builder.command()) : null
            }
            this.setStatus(TaskStatus.COMPLETED)
            throw new ProcessFailedException("Error submitting array job for execution", e)
        }
    }

    protected void setJobId(String jobId) {
        this.jobId = jobId
        array.eachWithIndex { handler, i ->
            ((GridTaskHandler)handler).setJobId(executor.getArrayTaskId(jobId, i))
        }
    }

    protected ProcessBuilder createProcessBuilder() {

        // -- log the submit command
        final cli = executor.getArraySubmitCommandLine()
        log.trace "submit array job > cli: ${cli}"

        // -- launch array job script
        new ProcessBuilder()
            .command( cli as String[] )
            .redirectErrorStream(true)
    }

    protected String processStart(ProcessBuilder builder, String arrayScript) {
        final process = builder.start()

        try {
            // -- pipe the array job script to the command stdin 
            log.trace "[${executor.name.toUpperCase()}] Submit array job >\n${arrayScript.indent()}"
            process.out << arrayScript
            process.out.close()

            // -- wait for the submission to complete
            final result = process.text
            final exitStatus = process.waitFor()
            final cmd = launchCmd0(builder, arrayScript)

            if( exitStatus )
                throw new ProcessNonZeroExitStatusException("Failed to submit array job to grid scheduler for execution", result, exitStatus, cmd)

            // -- return the process stdout
            return result
        }
        finally {
            // make sure to release all resources
            process.in.closeQuietly()
            process.out.closeQuietly()
            process.err.closeQuietly()
            process.destroy()
        }
    }

    protected String launchCmd0(ProcessBuilder builder, String arrayScript) {
        final cmd = CmdLineHelper.toLine(builder.command())

        new StringBuilder()
            .append("cat << 'LAUNCH_COMMAND_EOF' | ${cmd}\n")
            .append(arrayScript.trim())
            .append('\nLAUNCH_COMMAND_EOF\n')
            .toString()
    }

    @Override
    void kill() {
        executor.killTask(jobId)
    }

    @Override
    protected StringBuilder toStringBuilder(StringBuilder builder) {
        builder << "\n    array jobId: $jobId; "

        return super.toStringBuilder(builder)
    }

}
