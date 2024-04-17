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

import java.nio.file.LinkOption
import java.nio.file.Path
import java.nio.file.NoSuchFileException

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.exception.IllegalArityException
import nextflow.exception.MissingFileException
import nextflow.file.FileHelper
import nextflow.file.FilePatternSplitter
import nextflow.script.ProcessFileOutput
/**
 * Implements the collection of files from the work directory
 * of a task execution.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class TaskFileCollecter {

    private ProcessFileOutput param

    private TaskRun task

    private Path workDir

    TaskFileCollecter(ProcessFileOutput param, TaskRun task) {
        this.param = param
        this.task = task
        this.workDir = task.getTargetDir()
    }

    Object collect() {
        final List<Path> allFiles = []
        final filePatterns = param.getFilePatterns(task.context, workDir)
        boolean inputsExcluded = false

        for( String filePattern : filePatterns ) {
            List<Path> result = null

            final splitter = param.glob ? FilePatternSplitter.glob().parse(filePattern) : null
            if( splitter?.isPattern() ) {
                result = fetchResultFiles(filePattern, workDir)
                if( result && !param.includeInputs ) {
                    result = excludeStagedInputs(task, result)
                    log.trace "Process ${task.lazyName()} > after removing staged inputs: ${result}"
                    inputsExcluded |= (result.size()==0)
                }
            }
            else {
                final path = param.glob ? splitter.strip(filePattern) : filePattern
                final file = workDir.resolve(path)
                final exists = checkFileExists(file)
                if( exists )
                    result = List.of(file)
                else
                    log.debug "Process `${task.lazyName()}` is unable to find [${file.class.simpleName}]: `$file` (pattern: `$filePattern`)"
            }

            if( result )
                allFiles.addAll(result)

            else if( !param.optional && (!param.arity || param.arity.min > 0) ) {
                def msg = "Missing output file(s) `$filePattern` expected by process `${task.lazyName()}`"
                if( inputsExcluded )
                    msg += " (note: input files are not included in the default matching set)"
                throw new MissingFileException(msg)
            }
        }

        if( !param.isValidArity(allFiles.size()) )
            throw new IllegalArityException("Incorrect number of output files for process `${task.lazyName()}` -- expected ${param.arity}, found ${allFiles.size()}")

        return allFiles.size()==1 && param.isSingle() ? allFiles[0] : allFiles
    }

    /**
     * Collect the file(s) matching the specified name or glob pattern
     * in the given task work directory.
     *
     * @param pattern
     * @param workDir
     */
    protected List<Path> fetchResultFiles(String pattern, Path workDir) {
        final opts = [
            relative: false,
            hidden: param.hidden ?: pattern.startsWith('.'),
            followLinks: param.followLinks,
            maxDepth: param.maxDepth,
            type: param.type ? param.type : ( pattern.contains('**') ? 'file' : 'any' )
        ]

        List<Path> files = []
        try {
            FileHelper.visitFiles(opts, workDir, pattern) { Path it -> files.add(it) }
        }
        catch( NoSuchFileException e ) {
            throw new MissingFileException("Cannot access directory: '$workDir'", e)
        }

        return files.sort()
    }

    /**
     * Remove each path in the given list whose name matches the name of
     * an input file for the specified {@code TaskRun}
     *
     * @param task
     * @param collectedFiles
     */
    protected List<Path> excludeStagedInputs(TaskRun task, List<Path> collectedFiles) {

        final List<String> allStagedFiles = task.inputFiles.collect { it.stageName }
        final List<Path> result = new ArrayList<>(collectedFiles.size())

        for( int i = 0; i < collectedFiles.size(); i++ ) {
            final file = collectedFiles.get(i)
            final relativeName = workDir.relativize(file).toString()
            if( !allStagedFiles.contains(relativeName) )
                result.add(file)
        }

        return result
    }

    protected boolean checkFileExists(Path file) {
        param.followLinks ? file.exists() : file.exists(LinkOption.NOFOLLOW_LINKS)
    }
}
