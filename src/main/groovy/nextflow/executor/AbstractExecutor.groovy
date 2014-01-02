package nextflow.executor

import java.nio.file.Path
import java.util.concurrent.CountDownLatch

import groovy.io.FileType
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.exception.MissingFileException
import nextflow.processor.FileHolder
import nextflow.processor.FileOutParam
import nextflow.processor.InParam
import nextflow.processor.TaskConfig
import nextflow.processor.TaskHandler
import nextflow.processor.TaskMonitor
import nextflow.processor.TaskRun

/**
 * Declares methods have to be implemented by a generic
 * execution strategy
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
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
        if( !monitor ) {
            monitor = session.dispatcher.getOrCreateMonitor(this.class) { createTaskMonitor() }
        }
    }

    /**
     * @return Create a new instance of the {@code TaskQueueHolder} component
     */
    abstract protected TaskMonitor createTaskMonitor()

    /**
     * @return A reference to the current {@code #queueHolder} object
     */
    TaskMonitor getTaskMonitor()  { monitor }

    /**
     * @return Create a new {@code TaskHandler} to manage the scheduling
     * actions for this task
     */
    abstract TaskHandler createTaskHandler(TaskRun task)

    /**
     * Submit a task for execution
     *
     * @param task
     * @return
     */
    TaskHandler submitTask( TaskRun task, boolean blocking ) {
        def handler = createTaskHandler(task)
        monitor.put(handler)
        try {
            // set a count down latch if the execution is blocking
            if( blocking )
                handler.latch = new CountDownLatch(1)
            // now submit the task for execution
            handler.submit()
            return handler
        }
        catch( Exception e ) {
            monitor.remove(handler)
            return null
        }
    }


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
                throw new MissingFileException("Missing output file: '$fileName' expected by task: ${taskName}")
            }
            return result
        }

        // scan to find the file with that name
        List files = []
        workDirectory.eachFileMatch(FileType.ANY, ~/$filePattern/ ) { files << it }

        if( !files ) {
            throw new MissingFileException("Missing output file(s): '$fileName' expected by task: ${taskName}")
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
        def result = new StringBuilder()
        def fileOutNames = []
        task.getOutputsByType(FileOutParam).keySet().each { FileOutParam param ->
            fileOutNames.addAll( param.name.split( param.separatorChar ) as List )
        }

        // create a bash script that will copy the out file to the working directory
        log.debug "Unstaging file names: $fileOutNames"
        if( fileOutNames ) {
            result << 'for item in "' << fileOutNames.unique().join(' ') << '"; do' << separatorChar
            result << 'for name in `ls $item 2>/dev/null`; do' << separatorChar
            result << '  cp $name ' << task.workDirectory.toString() << separatorChar
            result << 'done' << separatorChar
            result << 'done' << separatorChar
        }

        return result.toString()
    }

}
