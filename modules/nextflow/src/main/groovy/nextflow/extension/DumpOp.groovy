/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

import org.codehaus.groovy.runtime.InvokerHelper
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
            log.info "[$marker] " + ( renderer ? renderer.call(it) : InvokerHelper.inspect(it) )
            target.bind(it)
        }

        final done = { DataflowExtensions.close(target) }

        DataflowHelper.subscribeImpl(source, [onNext: next, onComplete: done])
        return target
    }
}
