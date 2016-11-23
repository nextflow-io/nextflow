package nextflow.container
import java.nio.file.Path

import groovy.transform.PackageScope
import nextflow.util.Escape
import nextflow.util.MemoryUnit
import nextflow.util.PathTrie
/**
 * Base class for container wrapper builders
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
abstract class ContainerBuilder {

    final protected List env = []

    final protected List<Path> mounts = []

    protected List<String> runOptions = []

    protected List<String> engineOptions = []

    protected String cpus

    protected String memory

    protected String temp

    protected String image


    ContainerBuilder addRunOptions(String str) {
        runOptions.add(str)
        return this
    }

    ContainerBuilder addEngineOptions(String str) {
        engineOptions.add(str)
        return this
    }

    ContainerBuilder setCpus( String value ) {
        this.cpus = value
        return this
    }

    ContainerBuilder setMemory( value ) {
        if( value instanceof MemoryUnit )
            this.memory = "${value.toMega()}m"

        else if( value instanceof String )
            this.memory = value

        else
            throw new IllegalArgumentException("Not a supported memory value")

        return this
    }

    ContainerBuilder setName(String name) {

    }

    ContainerBuilder setTemp( String value ) {
        this.temp = value
        return this
    }

    abstract ContainerBuilder params( Map config )

    abstract ContainerBuilder build(StringBuilder result)

    abstract String getRunCommand()

    String getKillCommand() { return null }

    String getRemoveCommand() { return null }

    StringBuilder appendHelpers( StringBuilder wrapper ) {
        return wrapper
    }

    StringBuilder appendRunCommand( StringBuilder wrapper ) {
        wrapper << runCommand
    }

    ContainerBuilder build() {
        build(new StringBuilder())
    }

    ContainerBuilder addEnv( entry ) {
        env.add(entry)
        return this
    }

    ContainerBuilder addMount( Path path ) {
        if( path )
            mounts.add(path)
        return this
    }

    ContainerBuilder addMountForInputs( Map<String,Path> inputFiles ) {

        mounts.addAll( inputFilesToPaths(inputFiles) )
        return this
    }

    @PackageScope
    static List<Path> inputFilesToPaths( Map<String,Path> inputFiles ) {

        def List<Path> files = []
        inputFiles.each { name, storePath ->

            def path = storePath.getParent()
            if( path ) files << path

        }
        return files

    }

    /**
     * Given a normalised shell script (starting with a she-bang line)
     * replace the first token on the first line with a docker run command
     *
     * @param script
     * @param docker
     * @return
     */
    String addContainerRunCommand( ContainerScriptTokens script ) {

        final result = new ArrayList<String>(script.lines)
        final i = script.index
        final main = result[i].trim()
        final p = main.indexOf(' ')
        result[i] = ( p != -1
                ? runCommand + main.substring(p)
                : runCommand )
        result.add('')
        return result.join('\n')
    }

    /**
     * Get the env command line option for the give environment definition
     *
     * @param env
     * @param result
     * @return
     */
    protected CharSequence makeEnv( env, StringBuilder result = new StringBuilder() ) {
        // append the environment configuration
        if( env instanceof File ) {
            env = env.toPath()
        }
        if( env instanceof Path ) {
            result << '-e "BASH_ENV=' << Escape.path(env) << '"'
        }
        else if( env instanceof Map ) {
            short index = 0
            for( Map.Entry entry : env.entrySet() ) {
                if( index++ ) result << ' '
                result << ("-e \"${entry.key}=${entry.value}\"")
            }
        }
        else if( env instanceof String && env.contains('=') ) {
            result << '-e "' << env << '"'
        }
        else if( env ) {
            throw new IllegalArgumentException("Not a valid environment value: $env [${env.class.name}]")
        }

        return result
    }


    /**
     * Get the volumes command line options for the given list of input files
     *
     * @param mountPaths
     * @param binDir
     * @param result
     * @return
     */
    protected CharSequence makeVolumes(List<Path> mountPaths, StringBuilder result = new StringBuilder() ) {

        // find the longest commons paths and mount only them
        def trie = new PathTrie()
        mountPaths.each { trie.add(it) }

        def paths = trie.longest()
        paths.each{ if(it) result << "-v ${Escape.path(it)}:${Escape.path(it)} " }

        // -- append by default the current path -- this is needed when `scratch` is set to true
        result << '-v "$PWD":"$PWD"'

        return result
    }

}
