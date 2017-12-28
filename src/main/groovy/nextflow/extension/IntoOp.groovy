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

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.expression.DataflowExpression
import groovyx.gpars.dataflow.operator.ChainWithClosure
import groovyx.gpars.dataflow.operator.CopyChannelsClosure
import groovyx.gpars.dataflow.operator.DataflowEventAdapter
import groovyx.gpars.dataflow.operator.DataflowProcessor
import nextflow.Channel
import nextflow.Global
import nextflow.Session
/**
 * Implements the {@link DataflowExtensions#into} operators logic
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class IntoOp {

    private DataflowReadChannel source

    private List<DataflowWriteChannel> outputs

    private Session session = (Session)Global.session


    IntoOp( DataflowReadChannel source, List<DataflowWriteChannel> targets ) {
        assert source
        assert targets

        this.source = source
        this.outputs = targets
    }

    IntoOp( DataflowReadChannel source, int n ) {
        assert source
        assert n

        def targets = new ArrayList(n)
        for( int i=0; i<n; i++ )
            targets << new DataflowQueue()

        this.source = source
        this.outputs = targets
    }

    IntoOp( DataflowReadChannel source, Closure holder ) {
        assert source
        assert holder

        final names = CaptureProperties.capture(holder)
        final binding = Global.session.binding
        if( !names )
            throw new IllegalArgumentException("Missing target channel names in `into` operator")
        if( names.size() == 1 )
            log.warn("The `into` operator should be used to connect two or more target channels -- consider to replace it with `.set { ${names[0]} }`")

        def targets = []
        names.each { identifier ->
            def channel = DataflowExtensions.newChannelBy(source)
            targets.add(channel)
            binding.setVariable(identifier, channel)
        }

        this.source = source
        this.outputs = targets
    }

    List<DataflowWriteChannel> getOutputs() { outputs }

    IntoOp apply() {

        final params = [:]
        params.inputs = [source]
        params.outputs = outputs
        params.listeners = createListener()

        DataflowHelper.newOperator(params, new ChainWithClosure(new CopyChannelsClosure()))

        return this
    }

    private createListener() {

        final stopOnFirst = source instanceof DataflowExpression
        final listener = new DataflowEventAdapter() {
            @Override
            void afterRun(DataflowProcessor processor, List<Object> messages) {
                if( !stopOnFirst ) return
                // -- terminate the process
                processor.terminate()
                // -- close the output channels
                for( def it : outputs ) {
                    if( !(it instanceof DataflowExpression))
                        it.bind(Channel.STOP)

                    else if( !(it as DataflowExpression).isBound() )
                        it.bind(Channel.STOP)

                }
            }

            @Override
            public boolean onException(final DataflowProcessor processor, final Throwable e) {
                log.error("@unknown", e)
                session.abort(e)
                return true;
            }
        }

        return [listener]
    }

}
