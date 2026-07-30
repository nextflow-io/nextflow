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

package nextflow.extension

import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.Session

/**
 * Interface for nf-lineage extensions.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io
 */
interface LinExtension {

    static final Map PARAMS = [
        label: [List,String,GString],
        taskRun: [String,GString],
        workflowRun: [String,GString],
    ]

    /**
     * Query Lineage metadata to get files produced by tasks, workflows or annotations.
     *
     * @param session Nextflow Session
     * @param channel Channel to publish the Lineage Ids matching the query params
     * @param params Parameters for the lineage metadata query
     */
    abstract void fromLineage(Session session, DataflowWriteChannel channel, Map<String,?> params)
}
