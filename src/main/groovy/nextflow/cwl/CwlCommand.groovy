package nextflow.cwl

/**
 * Wrap a CWL CommandLineTool command line
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CwlCommand extends Closure {

    String command

    CwlCommand(String command) {
        super(null)
        this.command = command
    }

    @Override
    Object call() {
        return command
    }
}
