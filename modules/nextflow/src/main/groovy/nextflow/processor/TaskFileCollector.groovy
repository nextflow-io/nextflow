/*
 * Copyright 2013-2025, Seqera Labs
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
import nextflow.exception.MissingFileException
import nextflow.file.FileHelper
import nextflow.file.FilePatternSplitter
/**
 * Implements the collection of output files from a task
 * work directory.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class TaskFileCollector {

    private static final Map<String,?> DEFAULT_OPTS = [followLinks: true, glob: true]

    private List<String> filePatterns

    private Map<String,?> opts

    private TaskRun task

    private Path workDir

    TaskFileCollector(List<String> filePatterns, Map<String,?> opts, TaskRun task) {
        this.filePatterns = filePatterns
        this.opts = DEFAULT_OPTS + opts
        this.task = task
        this.workDir = task.getTargetDir()
    }

    List<Path> collect() {
        final List<Path> allFiles = []
        boolean inputsExcluded = false

        for( String filePattern : filePatterns ) {
            List<Path> result = null

            final splitter = opts.glob ? FilePatternSplitter.glob().parse(filePattern) : null
            if( splitter?.isPattern() ) {
                result = fetchResultFiles(filePattern, workDir)
                if( result && !opts.includeInputs ) {
                    result = excludeStagedInputs(result)
                    log.trace "Process ${task.lazyName()} > after removing staged inputs: ${result}"
                    inputsExcluded |= (result.size()==0)
                }
            }
            else {
                final path = opts.glob ? splitter.strip(filePattern) : filePattern
                final file = workDir.resolve(path)
                final exists = checkFileExists(file)
                if( exists )
                    result = List.of(file)
                else
                    log.debug "Process `${task.lazyName()}` is unable to find [${file.class.simpleName}]: `$file` (pattern: `$filePattern`)"
            }

            if( result ) {
                allFiles.addAll(result)
            }
            else if( !opts.optional ) {
                def msg = "Missing output file(s) `$filePattern` expected by process `${task.lazyName()}`"
                if( inputsExcluded )
                    msg += " (note: input files are not included in the default matching set)"
                throw new MissingFileException(msg)
            }
        }

        return allFiles
    }

    /**
     * Collect the file(s) matching the specified name or glob pattern
     * in the given task work directory.
     *
     * @param pattern
     * @param workDir
     */
    protected List<Path> fetchResultFiles(String pattern, Path workDir) {
        final opts = visitOptions(pattern)

        List<Path> files = []
        try {
            FileHelper.visitFiles(opts, workDir, pattern) { Path it -> files.add(it) }
        }
        catch( NoSuchFileException e ) {
            throw new MissingFileException("Cannot access directory: '$workDir'", e)
        }

        return files.sort()
    }

    protected Map<String,?> visitOptions(String pattern) {
        return [
            relative: false,
            hidden: opts.hidden ?: pattern.startsWith('.'),
            followLinks: opts.followLinks,
            maxDepth: opts.maxDepth,
            type: opts.type ? opts.type : ( pattern.contains('**') ? 'file' : 'any' )
        ]
    }

    /**
     * Remove each path in the given list whose name matches the name of
     * an input file for the specified {@code TaskRun}
     *
     * @param collectedFiles
     */
    protected List<Path> excludeStagedInputs(List<Path> collectedFiles) {

        final List<String> allStagedFiles = task.getStagedInputs()
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
        return opts.followLinks ? file.exists() : file.exists(LinkOption.NOFOLLOW_LINKS)
    }
}
