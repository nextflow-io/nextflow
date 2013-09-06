package nextflow.executor
import groovy.io.FileType
import groovy.util.logging.Slf4j
import nextflow.exception.MissingFileException
import nextflow.processor.EnvInParam
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
abstract class AbstractExecutor {

    /*
     * The object holding the configuration declared by this task
     */
    public TaskConfig taskConfig

    /**
     * Execute the specified task shell script
     *
     * @param script The script string to be execute, e.g. a BASH script
     * @return {@code TaskDef}
     */
    abstract void launchTask( TaskRun task )

    /**
     * The file which contains the stdout produced by the executed task script
     *
     * @param task The user task to be executed
     * @return The absolute file to the produced script output
     */
    abstract getStdOutFile( TaskRun task )


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

        // the '-' stands for the script stdout, save to a file
        if( fileName == '-' ) {
            return getStdOutFile(task)
        }

        // replace any wildcards characters
        // TODO give a try to http://code.google.com/p/wildcard/  -or- http://commons.apache.org/io/
        String filePattern = fileName.replace("?", ".?").replace("*", ".*")

        // when there's not change in the pattern, try to find a single file
        if( filePattern == fileName ) {
            def result = new File(task.workDirectory,fileName)
            if( !result.exists() ) {
                throw new MissingFileException("Missing output file: '$fileName' expected by task: ${task.name}")
            }
            return result
        }

        // scan to find the file with that name
        List files = []
        task.workDirectory.eachFileMatch(FileType.FILES, ~/$filePattern/ ) { File it -> files << it}
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
    def void createEnvironmentFile( TaskRun task, File target ) {
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



}
