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
import nextflow.extension.LinExtension
import nextflow.lineage.serde.LinSerializable

import static nextflow.lineage.fs.LinPath.*

/**
 * Lineage channel extensions
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@CompileStatic
@Slf4j
class LinExtensionImpl implements LinExtension {

    @Override
    Object lineage(Session session, String lid) {
        final store = getStore(session)
        return LinUtils.getMetadataObject(store, new URI(lid))
    }

    @Override
    void queryLineage(Session session, DataflowWriteChannel channel, Map<String, String> params) {
        new LinPropertyValidator().validateQueryParams(params)
        final store = getStore(session)
        emitSearchResults(channel, store.search(params))
        channel.bind(Channel.STOP)
    }

    protected LinStore getStore(Session session) {
        final store = LinStoreFactory.getOrCreate(session)
        if( !store ) {
            throw new Exception("Lineage store not found - Check Nextflow configuration")
        }
        return store
    }

    private void emitSearchResults(DataflowWriteChannel channel, Map<String, LinSerializable> results) {
        if( !results ) {
            return
        }
        results.keySet().forEach { channel.bind(asUriString(it)) }
    }
}
