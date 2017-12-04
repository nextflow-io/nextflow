package nextflow.cwl
import java.nio.file.Path

import nextflow.file.FileHelper
import nextflow.script.BaseScript
/**
 * CWL task entry point.
 *
 * @See {@link BaseScript#cwlTask(java.lang.Object)}
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CwlTask extends Closure {

    private BaseScript owner

    private Path path

    CwlTask( Object owner, path ) {
        super(owner)
        this.owner = (BaseScript)owner
        this.path = path instanceof URI ? FileHelper.asPath((URI)path) : FileHelper.asPath(path.toString())
    }

    @Override
    Object call() {
        def cwl = new CwlRunner(owner)
        cwl.run(path)
    }

    @Override
    Object call(Object... args) {
        throw new UnsupportedOperationException()
    }

    @Override
    Object call(Object args) {
        if( !(args instanceof Map) ) {
            throw new IllegalArgumentException("Not a valid cwlTask argument type: $args")
        }

        def cwl = new CwlRunner(owner, args)
        cwl.run(path)
    }

}
