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
import nextflow.lineage.fs.LinPathFactory
import nextflow.lineage.model.v1beta1.FileOutput

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
    void fromLineage(Session session, DataflowWriteChannel channel, Map<String,?> opts) {
        final queryParams = buildQueryParams(opts)
        log.trace("Querying lineage with params: $queryParams")
        new LinPropertyValidator().validateQueryParams(queryParams.keySet())
        final store = getStore(session)
        try( def stream = store.search(queryParams) ){
            stream.forEach { channel.bind( LinPathFactory.create( asUriString(it) ) ) }
        }
        channel.bind(Channel.STOP)
    }

    private static Map<String, List<String>> buildQueryParams(Map<String,?> opts) {
        final queryParams = [type: [FileOutput.class.simpleName] ]
        if( opts.workflowRun )
            queryParams['workflowRun'] = [opts.workflowRun as String]
        if( opts.taskRun )
            queryParams['taskRun'] = [opts.taskRun as String]
        if( opts.label ) {
            if( opts.label instanceof List )
                queryParams['labels'] = opts.label as List<String>
            else
                queryParams['labels'] = [ opts.label.toString() ]
        }
        return queryParams
    }

    protected LinStore getStore(Session session) {
        final store = LinStoreFactory.getOrCreate(session)
        if( !store ) {
            throw new Exception("Lineage store not found - Check Nextflow configuration")
        }
        return store
    }
}
