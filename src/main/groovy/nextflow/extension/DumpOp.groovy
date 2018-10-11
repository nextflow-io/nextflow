/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
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
