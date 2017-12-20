package nextflow.extension
import static nextflow.util.CheckHelper.checkParams

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowReadChannel
import nextflow.Global
import nextflow.Session

/**
 * Implements channel `dump` operator. It prints the content of a channel
 * only when the `-dump-channels` command line option is specified otherwise
 * it is ignored.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class DumpOp {

    static final private Map PARAMS_DUMP = [tag: String]

    private DataflowReadChannel source

    private Closure<String> renderer

    private List<String> dumpNames

    protected String tag

    DumpOp(DataflowReadChannel source, Map opts, Closure<String> renderer = null) {
        checkParams('dump', opts, PARAMS_DUMP)
        this.source = source
        this.tag = opts.tag
        this.renderer = renderer
        this.dumpNames = (Global.session as Session).getDumpChannels()
    }

    /** Only for testing -- do not use */
    protected DumpOp() {}


    boolean isEnabled() {
        if( !dumpNames )
            return false

        final matcher = tag ?: ''
        dumpNames
            .collect { it.replace('*','.*') }
            .find { matcher ==~ /$it/}
    }

    DataflowReadChannel apply() {

        if( !isEnabled() ) {
            return source
        }

        final target = DataflowExtensions.newChannelBy(source)

        final next = {
            def marker = 'DUMP'
            if( tag ) marker += ": $tag"
            log.info "[$marker] " + ( renderer ? renderer.call(it) : String.valueOf(it) )
            target.bind(it)
        }

        final done = { DataflowExtensions.close(target) }

        DataflowHelper.subscribeImpl(source, [onNext: next, onComplete: done])
        return target
    }
}
