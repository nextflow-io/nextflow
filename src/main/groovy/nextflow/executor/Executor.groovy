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
import nextflow.processor.TaskHandler
import nextflow.processor.TaskMonitor
import nextflow.processor.TaskRun
import nextflow.script.FileOutParam
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


}
