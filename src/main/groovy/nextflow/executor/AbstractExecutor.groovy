package nextflow.executor

import java.nio.file.Path

import groovy.io.FileType
import groovy.util.logging.Slf4j
import nextflow.exception.MissingFileException
import nextflow.processor.EnvInParam
import nextflow.processor.FileInParam
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
        task.workDirectory.eachFileMatch(FileType.FILES, ~/$filePattern/ ) { Path it -> files << it}
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
    def void createEnvironmentFile( TaskRun task, Path target ) {
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
    def String stagingFilesScript( Map<FileInParam,Object> inputs ) {
        assert inputs != null

        def count = 0
        def delete = []
        def links = []
        inputs.each { param, obj ->

            def files = (obj instanceof Collection ? obj : [obj]).collect {
                normalizeInputToFile(it, "input.${++count}")
            }

            def names = expandWildcards(param.name, files)

            // delete all previous files with the same name
            names.each {
                delete << "rm -f ${it}"
            }

            // link them
            names.eachWithIndex { String entry, int i ->
                links << stageInputFileScript(files[i], entry)
            }
        }
        links << '' // just to have new-line at the end of the script

        // return a big string containing the command
        return (delete + links).join('; ')
    }

    /**
     * An input file parameter can be provided with any value other than a file.
     * This function normalize a generic value to a {@code Path} create a temporary file
     * in the for it.
     *
     * @param input The input value
     * @param altName The name to be used when a temporary file is created.
     * @return The {@code Path} that will be staged in the task working folder
     */
    protected Path normalizeInputToFile( Object input, String altName ) {

        if( input instanceof Path ) {
            return input
        }

        def result = taskConfig.getOwnerScript().tempFile(altName)
        result.text = input?.toString() ?: ''
        return result
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
    protected String stageInputFileScript( Path path, String targetName ) {
        "ln -s ${path.toAbsolutePath()} $targetName"
    }

    /**
     * An input file name may contain wildcards characters which have to be handled coherently
     * given the number of files specified.
     *
     * @param name A file name with may contain a wildcard character star {@code *} or question mark {@code ?}.
     *  Only one occurrence can be specified for star or question mark widlcards.
     *
     * @param value Any value that have to be managed as an input files. Values other than {@code Path} are converted
     * to a string value, using the {@code #toString} method and saved in the local file-system. Value of type {@code Collection}
     * are expanded to multiple values accordingly.
     *
     * @return
     */
    protected List<String> expandWildcards( String name, Object value ) {
        assert name
        assert value != null

        def result = []
        if( name == '*' ) {
            def files = value instanceof Collection ? value : [value]
            files.each {
                if( it instanceof Path ) { result << it.name }
                else throw new IllegalArgumentException("Not a valid value argument for 'expandWildcards' method: $it")
            }
            return result
        }

        // no wildcards in the file name
        else if( !name.contains('*') && !name.contains('?') ) {

            /*
             * The name do not contain any wildcards *BUT* when multiple files are provide
             * it is managed like having a 'start' at the end of the file name
             */
            if( value instanceof Collection ) {
                name += '*'
            }
            else {
                // just return that name
                return [name]
            }
        }

        /*
         * The star wildcard: when a single item is provided, it is simply ignored
         * When a collection of files is provided, the name is expanded to the index number
         */
        if( name.contains('*') ) {
            if( value instanceof Collection && value.size()>1 ) {
                def count = 1
                value.each {
                    result << name.replace('*', (count++).toString())
                }
            }
            else {
                // there's just one value, remove the 'star' wildcards
                result << name.replace('*','')
            }
        }

        /*
         * The question mark wildcards *always* expand to an index number
         * as long as are the number of question mark characters
         */
        else if( name.contains('?') ) {
            def files = value instanceof Collection ? value : [value]
            def count = 1
            files.each {
                String match = (name =~ /\?+/)[0]
                def replace = (count++).toString().padLeft(match.size(), '0')
                def fileName = name.replace(match, replace)
                result << fileName
            }

        }

        // not a valid condition
        else {
            throw new IllegalStateException("Invalid file expansion for name: '$name'")
        }

        return result
    }





}
