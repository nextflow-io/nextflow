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
 */
package nextflow.processor

import java.lang.reflect.InvocationTargetException
import java.nio.file.FileSystems
import java.nio.file.NoSuchFileException
import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import nextflow.exception.FailedGuardException
import nextflow.exception.ProcessEvalException
import nextflow.exception.ShowOnlyExceptionMessage
import nextflow.plugin.Plugins
import nextflow.processor.tip.TaskTipProvider
import nextflow.util.LoggerHelper
/**
 * Implement formatting of standard task errors.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class TaskErrorFormatter {

    /**
     * Format the error message for an eval output error.
     *
     * @param message
     * @param error
     * @param task
     */
    public List<String> formatCommandError(List<String> message, ProcessEvalException error, TaskRun task) {
        // compose a readable error message
        message << formatErrorCause(error)

        // - print the executed command
        message << "Command executed:\n"
        for( final line : error.command.stripIndent(true)?.trim()?.readLines() )
            message << "  ${line}".toString()

        // - the exit status
        message << "\nCommand exit status:\n  ${error.status}".toString()

        // - the tail of the process stdout
        message << "\nCommand output:"

        final lines = error.output.readLines()
        if( lines.size() == 0 )
            message << "  (empty)"
        for( final line : lines )
            message << "  ${stripWorkDir(line, task.workDir)}".toString()

        if( task?.workDir )
            message << "\nWork dir:\n  ${task.workDirStr}".toString()

        return message
    }

    /**
     * Format the error message for a guard error (e.g. `when:`).
     *
     * @param message
     * @param error
     * @param task
     */
    public List<String> formatGuardError(List<String> message, FailedGuardException error, TaskRun task) {
        // compose a readable error message
        message << formatErrorCause(error)

        if( error.source ) {
            message << "\nWhen block:"
            for( final line : error.source.stripIndent(true).readLines() )
                message << "  ${line}".toString()
        }

        if( task?.workDir )
            message << "\nWork dir:\n  ${task.workDirStr}".toString()

        return message
    }

    /**
     * Format the error message for a task error.
     *
     * @param message
     * @param error
     * @param task
     */
    public List<String> formatTaskError(List<String> message, Throwable error, TaskRun task) {

        // prepend a readable error message
        message << formatErrorCause( error )

        // task with `script:` block
        if( task?.script ) {
            // -- print the executed command
            message << "Command executed${task.template ? " [$task.template]": ''}:\n".toString()
            for( final line : task.script?.stripIndent(true)?.trim()?.readLines() )
                message << "  ${line}".toString()

            // -- the exit status
            message << "\nCommand exit status:\n  ${task.exitStatus != Integer.MAX_VALUE ? task.exitStatus : '-'}".toString()

            // -- the tail of the process stdout
            message << "\nCommand output:"
            final max = 50
            final stdoutLines = task.dumpStdout(max)
            if( stdoutLines.size() == 0 )
                message << "  (empty)"
            for( final line : stdoutLines )
                message << "  ${stripWorkDir(line, task.workDir)}".toString()

            // -- the tail of the process stderr
            final stderrLines = task.dumpStderr(max)
            if( stderrLines ) {
                message << "\nCommand error:"
                for( final line : stderrLines )
                    message << "  ${stripWorkDir(line, task.workDir)}".toString()
            }

            // -- this is likely a task wrapper issue
            else if( task.exitStatus != 0 ) {
                final logLines = task.dumpLogFile(max)
                if( logLines ) {
                    message << "\nCommand log:"
                    for( final line : logLines )
                        message << "  ${stripWorkDir(line, task.workDir)}".toString()
                }
            }
        }

        // task with `exec:` block
        else if( task?.source ) {
            message << "Source block:"
            for( final line : task.source.stripIndent(true).readLines() )
                message << "  ${line}".toString()
        }

        // append work dir if present
        if( task?.workDir )
            message << "\nWork dir:\n  ${task.workDirStr}".toString()

        // append container image if present
        if( task?.isContainerEnabled() )
            message << "\nContainer:\n  ${task.container}".toString()

        // append tip
        message << suggestTip(message)

        return message
    }

    /**
     * Format the cause of an error.
     *
     * @param error
     */
    public String formatErrorCause(Throwable error) {

        final result = new StringBuilder()
        result.append('\nCaused by:\n')

        final message = error instanceof ShowOnlyExceptionMessage || !error.cause
            ? err0(error)
            : err0(error.cause)

        for( final line : message.readLines() )
            result.append('  ').append(line).append('\n')

        return result.append('\n').toString()
    }

    private static String err0(Throwable e) {
        final err = e instanceof InvocationTargetException ? e.targetException : e

        if( err instanceof NoSuchFileException ) {
            return "No such file or directory: $err.message"
        }

        if( err instanceof MissingPropertyException ) {
            def name = err.property ?: LoggerHelper.getDetailMessage(err)
            def result = "No such variable: ${name}"
            def details = LoggerHelper.findErrorLine(err)
            if( details )
                result += " -- Check script '${details[0]}' at line: ${details[1]}"
            return result
        }

        def result = err.message ?: err.toString()
        def details = LoggerHelper.findErrorLine(err)
        if( details ) {
            result += " -- Check script '${details[0]}' at line: ${details[1]}"
        }
        return result
    }

    private static String stripWorkDir(String line, Path workDir) {
        if( workDir == null )
            return line

        if( workDir.fileSystem != FileSystems.default )
            return line

        return workDir ? line.replace(workDir.toString() + '/', '') : line
    }

    private static String suggestTip(List<String> message) {
        try {
            return "\nTip: ${getTipProvider().suggestTip(message)}"
        }
        catch (Exception e) {
            log.debug "Unable to get tip for task message: $message", e
            return ''
        }
    }

    @Memoized
    private static TaskTipProvider getTipProvider() {
        final provider = Plugins.getPriorityExtensions(TaskTipProvider).find(it-> it.enabled())
        if( !provider )
            throw new IllegalStateException("Unable to find any tip provider")
        return provider
    }

}
