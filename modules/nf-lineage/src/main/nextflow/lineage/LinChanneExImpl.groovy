/*
 * Copyright 2013-2025, Seqera Labs
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

package nextflow.lineage

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.Channel
import nextflow.Session
import nextflow.extension.LinChannelEx
import nextflow.lineage.fs.LinPath
import nextflow.lineage.fs.LinPathFactory
import nextflow.lineage.serde.LinSerializable

/**
 * Lineage channel extensions
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@CompileStatic
@Slf4j
class LinChanneExImpl implements LinChannelEx{

    void viewLineage(Session session, DataflowWriteChannel channel, URI uri) {
        final store = getStore(session)
        emitResults(channel, LinUtils.query(store, uri))
        channel.bind(Channel.STOP)
    }

    void queryLineage(Session session, DataflowWriteChannel channel, String query) {
        final store = getStore(session)
        emitSearchResults(channel, store.search(query))
        channel.bind(Channel.STOP)
    }


    protected LinStore getStore(Session session){
        final store = LinStoreFactory.getOrCreate(session)
        if( !store ) {
            throw new Exception("Lineage store not found - Check Nextflow configuration")
        }
        return store
    }

    private static void emitResults(DataflowWriteChannel channel, Collection results){
        if( !results ) {
            return
        }
        // Remove nested collections of a single element
        if( results.size() == 1 ) {
            final entry = results[0]
            if( entry instanceof Collection ) {
                emitResults(channel, entry)
            } else {
                channel.bind(LinUtils.encodeSearchOutputs(entry))
            }
        } else
            results.forEach { channel.bind(LinUtils.encodeSearchOutputs(it)) }
    }

    private void emitSearchResults(DataflowWriteChannel channel, Map<String, LinSerializable> results) {
        if( !results ) {
            return
        }
        results.keySet().forEach { channel.bind(LinPathFactory.create(LinPath.LID_PROT + it)) }
    }
}
