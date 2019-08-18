/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package nextflow.extension

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.Global
import nextflow.Session
import org.codehaus.groovy.runtime.InvokerHelper
import static nextflow.util.CheckHelper.checkParams
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

    private Session session = (Global.session as Session)

    private DataflowReadChannel source

    private Closure<String> renderer

    private List<String> dumpNames

    protected String tag

    DumpOp(Map opts, Closure<String> renderer) {
        checkParams('dump', opts, PARAMS_DUMP)
        this.source = source
        this.tag = opts.tag
        this.renderer = renderer
        this.dumpNames = session.getDumpChannels()
    }

    DumpOp setSource( DataflowWriteChannel source ) {
        this.source = CH.getReadChannel(source)
        return this
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

    DataflowWriteChannel apply() {

        if( !isEnabled() ) {
            if( source instanceof DataflowWriteChannel )
                return (DataflowWriteChannel)source
            throw new IllegalArgumentException("Illegal dump operator source channel")
        }

        final target = CH.createBy(source)
        final events = new HashMap(2)
        events.onNext = {
            def marker = 'DUMP'
            if( tag ) marker += ": $tag"
            log.info "[$marker] " + ( renderer ? renderer.call(it) : InvokerHelper.inspect(it) )
            target.bind(it)
        }

        events.onComplete = { CH.close0(target) }

        DataflowHelper.subscribeImpl(source, events)
        return target
    }
}
