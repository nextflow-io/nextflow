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
package nextflow.lineage

/**
 * Interface to log workflow executions and their corresponding Lineage IDs
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
interface LinHistoryLog {
    /**
     * Write a workflow execution linage history log record.
     *
     * @param name Workflow execution name.
     * @param sessionId Workflow session ID.
     * @param runLid Workflow run ID.
     */
    void write(String name, UUID sessionId, String runLid)

    /**
     * Get the store records in the Lineage History Log.
     *
     * @return List of stored lineage history records.
     */
    List<LinHistoryRecord> getRecords()

}
