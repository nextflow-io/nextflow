package nextflow.container
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
trait ContainerBuilder {

    abstract ContainerBuilder params( Map config )

    abstract String build(StringBuilder result)

    abstract getRunCommand()

    String getKillCommand() { return null }

    String getRemoveCommand() { return null }

    StringBuilder appendHelpers( StringBuilder wrapper ) {
        return wrapper
    }

    StringBuilder appendRunCommand( StringBuilder wrapper ) {
        wrapper << runCommand
    }

    String build() {
        build(new StringBuilder())
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

}
