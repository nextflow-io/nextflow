/*
 * Copyright 2013-2026, Seqera Labs
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

package nextflow.script

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowVariable
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.dataflow.ChannelImpl
import nextflow.dataflow.ValueImpl
import nextflow.util.RecordMap
import org.codehaus.groovy.runtime.InvokerHelper
/**
 * Utility functions for converting between v1 and v2 dataflow types.
 *
 * When a script enables the `nextflow.preview.types` flag, its workflows
 * are "typed workflows", otherwise they are "legacy workflows":
 *
 * - Typed workflows use v2 dataflow types: ChannelImpl, ValueImpl
 * - Legacy workflows use v1 dataflow types: DataflowBroadcast, DataflowVariable
 *
 * While a given script must be either typed or legacy, typed workflows
 * can be composed with legacy workflows and vise versa. In order to support
 * this, the calling workflow must be able to convert between v1 and v2 dataflow
 * types based on whether the callee is typed or legacy.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class DataflowTypeHelper {

    /**
     * Normalize an array of arguments to the appropriate dataflow
     * type based on the target context.
     *
     * @param args
     * @param typingEnabled
     */
    static Object[] normalizeArray(Object args, boolean typingEnabled) {
        return InvokerHelper.asArray(args).stream()
            .map(arg -> normalize(arg, typingEnabled))
            .toArray()
    }

    /**
     * Normalize a source value to the appropriate dataflow type
     * based on the target context.
     *
     * If a v2 dataflow type (ChannelImpl or ValueImpl) is being
     * passed to a legacy process or workflow, normalize the source
     * value by unwrapping it into a DataflowWriteChannel.
     *
     * If a v1 dataflow type (DataflowWriteChannel) is being passed
     * to a typed process or workflow, normalize the source value by
     * wrapping it as a ChannelImpl or ValueImpl.
     *
     * @param source
     * @param typingEnabled
     */
    static Object normalize(Object source, boolean typingEnabled) {
        return typingEnabled
            ? normalizeV2(source)
            : normalizeV1(source)
    }

    /**
     * Normalize a source value by converting v1 dataflow types
     * (DataflowWriteChannel, ChannelOut) to v2 dataflow types
     * (ChannelImpl, ValueImpl).
     *
     * @param source
     */
    static Object normalizeV2(Object source) {
        if( source instanceof ChannelOut )
            return normalizeMultiChannelV2(source)
        if( source instanceof DataflowVariable )
            return new ValueImpl(source)
        if( source instanceof DataflowWriteChannel )
            return new ChannelImpl(source)
        return source
    }

    private static Object normalizeMultiChannelV2(ChannelOut source) {
        final names = source.getNames()
        if( source.size() == 1 )
            return normalizeV2(source[0])
        final result = new HashMap<String,Object>()
        for( int i = 0; i < source.size(); i++ )
            result.put(names[i], normalizeV2(source[i]))
        return new RecordMap(result)
    }

    /**
     * Normalize a source value by converting v2 dataflow types
     * (ChannelImpl, ValueImpl) to v1 dataflow types (DataflowWriteChannel,
     * ChannelOut).
     *
     * @param source
     */
    static Object normalizeV1(Object source) {
        if( source instanceof ChannelImpl )
            return source.getSource()
        if( source instanceof ValueImpl )
            return source.getSource()
        return source
    }

}
