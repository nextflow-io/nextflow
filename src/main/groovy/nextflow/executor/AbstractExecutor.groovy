/*
 * Copyright (c) 2013-2014, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2014, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
 */

package nextflow.executor
import java.nio.file.Path

import groovy.io.FileType
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.exception.MissingFileException
import nextflow.processor.FileHolder
import nextflow.processor.TaskConfig
import nextflow.processor.TaskHandler
import nextflow.processor.TaskMonitor
import nextflow.processor.TaskRun
import nextflow.script.InParam
import nextflow.script.ScriptType

/**
 * Declares methods have to be implemented by a generic
 * execution strategy
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@SupportedScriptTypes( [ScriptType.SCRIPTLET] )
abstract class AbstractExecutor {

    /**
     * The object holding the configuration declared by this task
     */
    TaskConfig taskConfig

    /**
     * The current session object
     */
    Session session

    /**
     * The executor simple name
     */
    String name

    /**
     * The queue holder that keep track of all tasks for this executor.
     */
    private TaskMonitor monitor


    /**
     * Let to post initialize the executor
     */
    def void init() {
        // -- skip if already assigned, this is only for testing purpose
        if( monitor ) return

        // -- get the reference to the monitor class for this executor
        monitor = session.dispatcher.getOrCreateMonitor(this.class) {
            log.info "[warm up] executor > $name"
            createTaskMonitor()
        }
    }


    /**
     * @return Create a new instance of the {@code TaskQueueHolder} component
     */
    abstract protected TaskMonitor createTaskMonitor()

    /**
     * @return A reference to the current {@code #queueHolder} object
     */
    @PackageScope
    TaskMonitor getTaskMonitor()  { monitor }

    /**
     * @return Create a new {@code TaskHandler} to manage the scheduling
     * actions for this task
     */
    abstract TaskHandler createTaskHandler(TaskRun task)


    /**
     * Collect the file(s) with the name specified, produced by the execution
     *
     * @param path The job working path
     * @param fileName The file name, it may include file name wildcards
     * @return The list of files matching the specified name
     */
    def collectResultFile( Path workDirectory, String fileName, String taskName ) {
        assert fileName
        assert workDirectory

        // replace any wildcards characters
        // TODO use newDirectoryStream here and eventually glob
        String filePattern = fileName.replace("?", ".?").replace("*", ".*")

        // when there's not change in the pattern, try to find a single file
        if( filePattern == fileName ) {
            def result = workDirectory.resolve(fileName)
            if( !result.exists() ) {
                throw new MissingFileException("Missing output file: '$fileName' expected by process: ${taskName}")
            }
            return result
        }

        // scan to find the file with that name
        List files = []
        try {
            workDirectory.eachFileMatch(FileType.ANY, ~/$filePattern/ ) { files << it }
        }
        catch( IOException e ) {
            throw new MissingFileException("Cannot access folder: '$workDirectory' expected by process: ${taskName}", e)
        }

        if( !files ) {
            throw new MissingFileException("Missing output file(s): '$fileName' expected by process: ${taskName}")
        }

        return files
    }



    /**
     * Given a map of the input file parameters with respective values,
     * create the BASH script to stage them into the task working space
     *
     * @param inputs An associative array mapping each {@code FileInParam} to the corresponding file (or generic value)
     * @return The BASH script to stage them
     */
    def String stagingFilesScript( Map<InParam, List<FileHolder>> inputs, String separatorChar = '\n') {
        assert inputs != null

        def delete = []
        def links = []
        inputs.each { param, files ->

            // delete all previous files with the same name
            files.each {
                delete << "rm -f ${it.stagePath.name}"
            }

            // link them
            files.each { FileHolder it ->
                links << stageInputFileScript( it.storePath, it.stagePath.name )
            }

        }
        links << '' // just to have new-line at the end of the script

        // return a big string containing the command
        return (delete + links).join(separatorChar)
    }

    /**
     * Stage the input file into the task working area. By default it creates a symlink
     * to the the specified path using {@code targetName} as name.
     * <p>
     *     An executor may override it to support a different staging mechanism
     *
     * @param path The {@code Path} to the file to be staged
     * @param targetName The name to be used in the task working directory
     * @return The script which will apply the staging for this file in the main script
     */
    String stageInputFileScript( Path path, String targetName ) {
        "ln -s ${path.toAbsolutePath()} $targetName"
    }

    /**
     * Creates the script to unstage the result output files from the scratch directory
     * to the shared working directory
     *
     * @param task The {@code TaskRun} executed
     * @param separatorChar The string to be used to separate multiple BASH statements (default: new line char)
     * @return The BASH script fragment to be used to copy the output files to the shared storage
     */
    String unstageOutputFilesScript( final TaskRun task, final String separatorChar = '\n' ) {

        // collect all the expected names (pattern) for files to be un-staged
        def result = []
        def fileOutNames = task.getOutputFilesNames()

        // create a bash script that will copy the out file to the working directory
        log.trace "Unstaging file names: $fileOutNames"
        if( fileOutNames ) {
            result << "mkdir -p ${task.getTargetDir()}"
            result << "for item in \"${fileOutNames.unique().join(' ')}\"; do"
            result << '  rsync -rRz $item ' + task.getTargetDir().toString()
            result << 'done'
        }

        return result.join(separatorChar)
    }

}
