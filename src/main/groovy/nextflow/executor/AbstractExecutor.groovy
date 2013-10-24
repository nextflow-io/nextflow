package nextflow.executor

import java.nio.file.Path

import groovy.io.FileType
import groovy.util.logging.Slf4j
import nextflow.exception.MissingFileException
import nextflow.processor.FileHolder
import nextflow.processor.FileInParam
import nextflow.processor.FileOutParam
import nextflow.processor.TaskConfig
import nextflow.processor.TaskHandler
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
    def collectResultFile( TaskRun task, String fileName ) {
        assert fileName
        assert task
        assert task.workDirectory

        // replace any wildcards characters
        // TODO use newDirectoryStream here and eventually glob
        String filePattern = fileName.replace("?", ".?").replace("*", ".*")

        // when there's not change in the pattern, try to find a single file
        if( filePattern == fileName ) {
            def result = task.workDirectory.resolve(fileName)
            if( !result.exists() ) {
                throw new MissingFileException("Missing output file: '$fileName' expected by task: ${task.name}")
            }
            return result
        }

        // scan to find the file with that name
        List files = []
        task.workDirectory.eachFileMatch(FileType.ANY, ~/$filePattern/ ) { files << it }

        if( !files ) {
            throw new MissingFileException("Missing output file(s): '$fileName' expected by task: ${task.name}")
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
    def String stagingFilesScript( Map<FileInParam, List<FileHolder>> inputs, String separatorChar = '\n') {
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
