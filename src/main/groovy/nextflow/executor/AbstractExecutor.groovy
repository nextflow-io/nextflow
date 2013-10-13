package nextflow.executor

import java.nio.file.Path

import groovy.io.FileType
import groovy.util.logging.Slf4j
import nextflow.exception.MissingFileException
import nextflow.processor.EnvInParam
import nextflow.processor.FileInParam
import nextflow.processor.FileHolder
import nextflow.processor.TaskConfig
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun

/**
 * Declares methods have to be implemented by a generic
 * execution strategy
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
abstract class AbstractExecutor<T extends ProcessHandler> {

    /**
     * The object holding the configuration declared by this task
     */
    TaskConfig taskConfig

    /**
     * Execute the specified task shell script
     *
     * @param script The script string to be execute, e.g. a BASH script
     * @return {@code TaskDef}
     */
    abstract void launchTask( TaskRun<T> task )

    /**
     * Check whenever the task has has started
     *
     * @param task
     * @return
     */
    abstract boolean checkStarted( TaskRun<T> task )

    /**
     * Check whenever the task has complete his job or it is still running
     *
     * @param task
     * @return
     */
    abstract boolean checkCompleted( TaskRun<T> task )

    /**
     * The file which contains the stdout produced by the executed task script
     *
     * @param task The user task to be executed
     * @return The absolute file to the produced script output
     */
    def getStdOutFile( TaskRun<T> task ) {
        task.handler.getOutputFile()
    }


    /**
     * Collect the file(s) with the name specified, produced by the execution
     *
     * @param path The job working path
     * @param fileName The file name, it may include file name wildcards
     * @return The list of files matching the specified name
     */
    def collectResultFile( TaskRun<T> task, String fileName ) {
        assert fileName
        assert task
        assert task.workDirectory

        // the '-' stands for the script stdout, save to a file
        if( fileName == '-' ) {
            return getStdOutFile(task)
        }

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
        task.workDirectory.eachFileMatch(FileType.FILES, ~/$filePattern/ ) { files << it }

        if( !files ) {
            throw new MissingFileException("Missing output file(s): '$fileName' expected by task: ${task.name}")
        }

        return files
    }

    /**
     * Save the environment into the specified file
     *
     * @param task The task instance which current environment needs to be stored
     * @param target The path to where save the task environment
     */
    def void createEnvironmentFile( TaskRun<T> task, Path target ) {
        assert task
        assert target

        // get the 'static' environment
        final environment = task.processor.getProcessEnvironment()

        // add the task input of type 'env'
        task.getInputsByType( EnvInParam ).each { param, value ->
            environment.put( param.name, value?.toString() )
        }

        // create the *bash* environment script
        target.text = TaskProcessor.bashEnvironmentScript(environment)

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


}
