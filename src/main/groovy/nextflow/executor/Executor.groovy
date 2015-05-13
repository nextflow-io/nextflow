/*
 * Copyright (c) 2013-2015, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2015, Paolo Di Tommaso and the respective authors.
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

import java.nio.file.NoSuchFileException
import java.nio.file.Path

import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.exception.MissingFileException
import nextflow.file.FileHelper
import nextflow.file.FileHolder
import nextflow.processor.TaskHandler
import nextflow.processor.TaskMonitor
import nextflow.processor.TaskRun
import nextflow.script.FileOutParam
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
abstract class Executor {

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

    private static final Map<Class,Boolean> registerFlag = [:]

    /**
     * Template method executed once for executor class
     */
    void register() { }

    /**
     * Allows to post-initialize the executor
     */
    void init() {
        log.debug "Initializing executor: $name"

        // -- skip if already assigned, this is only for testing purpose
        if( monitor )
            return

        // -- get the reference to the monitor class for this executor
        monitor = session.dispatcher.getOrCreateMonitor(this.class) {
            log.info "[warm up] executor > $name"
            createTaskMonitor()
        }

        // call the register template method
        if( !registerFlag[this.class] ) {
            log.debug "Invoke register for executor: ${name}"
            register()
            registerFlag[this.class] = true
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
     * @param workDir The job working path
     * @param fileName The file name, it may include file name wildcards
     * @return The list of files matching the specified name
     */
    def collectResultFile( Path workDir, String fileName, String taskName, FileOutParam param ) {
        assert fileName
        assert workDir

        List files = []
        def opts = collectResultOpts(param, fileName)
        // scan to find the file with that name
        try {
            FileHelper.visitFiles(opts, workDir, fileName) { Path it -> files.add(it) }
        }
        catch( NoSuchFileException e ) {
            throw new MissingFileException("Cannot access folder: '$workDir' expected by process: ${taskName}", e)
        }

        return files
    }

    protected Map collectResultOpts( FileOutParam param, String fileName ) {
        final opts = [:]
        opts.relative = false
        opts.hidden = param.hidden ?: fileName.startsWith('.')
        opts.followLinks = param.followLinks
        opts.maxDepth = param.maxDepth
        opts.type = param.type ? param.type : ( fileName.contains('**') ? 'file' : 'any' )
        return opts
    }


    def String stagingInputFilesScript( TaskRun task ) {
        stagingFilesScript(task.getInputFiles())
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
                delete << "rm -f '${it.stageName}'"
            }

            // link them
            files.each { FileHolder it ->
                links << stageInputFileScript( it.storePath, it.stageName )
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
        "ln -s '${path.toAbsolutePath()}' '$targetName'"
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
        def normalized = normalizeGlobStarPaths(fileOutNames)

        // create a bash script that will copy the out file to the working directory
        log.trace "Unstaging file path: $normalized"
        if( normalized ) {
            result << ""
            result << "mkdir -p ${task.getTargetDir()}"
            String strategy = task.getProcessor().getTaskConfig().unstageStrategy
            for( int i=0; i<normalized.size(); i++ ) {
                final path = normalized[i]
                final cmd = copyCommand(path, task.getTargetDir().toString(), strategy) + ' || true' // <-- add true to avoid it stops on errors
                result << cmd
            }
        }

        return result.join(separatorChar)
    }

    /**
     * Compose a `cp` or `mv` command given the source and target paths
     *
     * @param source
     * @param target
     * @param move When {@code true} use unix {@code mv} instead of {@code cp}
     * @return A shell copy or move command string
     */
    protected String copyCommand( String source, String target, String strategy = null) {

        def cmd
        if( strategy == 'copy' || !strategy )
            cmd = 'cp -fR'
        else if( strategy == 'move' )
            cmd = 'mv -f'
        else if( strategy == 'rsync' )
            return "rsync -rRl $source $target"
        else
            throw new IllegalArgumentException("Unknown un-stage strategy: $strategy")

        final p = source.lastIndexOf('/')
        if( p<=0 ) {
            return "$cmd $source $target"
        }

        def path  = new File(target,source.substring(0,p)).toString()
        return "mkdir -p $path && $cmd $source $path"
    }

    /**
     * Normalize path that contains glob star wildcards (i.e. double star character) since
     * they are not supported by plain BASH
     *
     * @param files
     * @return
     */
    protected List<String> normalizeGlobStarPaths( List<String> files ) {

        def result = []
        for( int i=0; i<files.size(); i++ ) {
            def item = removeGlobStar(files.get(i))
            if( !result.contains(item) )
                result.add(item)
        }

        return result
    }

    /**
     * Remove a glob star specifier returning the longest parent path
     *
     * <pre>
     * /some/data/&#42;&#42;/file.txt  ->   /some/data
     * /some/data/file.txt     ->   /some/data/file.txt
     * /some&#42;&#42;                 ->   *
     * </pre>
     *
     * @param path
     * @return
     */
    protected removeGlobStar(String path) {

        def p = path?.indexOf('**')
        if( p == -1 )
            return path

        def slash = -1
        for( int i=p-1; i>=0; i-- ) {
            if( path.charAt(i) == '/' as char) {
                slash = i
                break
            }
        }

        if( slash == -1 )
            return '*'

        path.substring(0,slash)
    }



}
