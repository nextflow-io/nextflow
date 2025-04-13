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
package nextflow.data.cid

/**
 * Interface to log workflow executions and their corresponding CIDs
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
interface CidHistoryLog {
    /**
     * Write a workflow execution CidHistoryLog record.
     *
     * @param name Workflow execution name.
     * @param sessionId Workflow session ID.
     * @param runCid Workflow run CID.
     * @param resultsCid Workflow results CID.
     */
    void write(String name, UUID sessionId, String runCid)

    /**
     * Updates the run CID for a given session ID.
     *
     * @param sessionId Workflow session ID.
     * @param runCid Workflow run CID.
     */
    void updateRunCid(UUID sessionId, String runCid)

    /**
     * Get the store records in the CidHistoryLog.
     *
     * @return List stored CIDHistoryRecords.
     */
    List<CidHistoryRecord> getRecords()

    /**
     * Get the record for a given
     * @param sessionId Workflow session ID.
     * @return CIDHistoryRecord for the given ID.
     */
    CidHistoryRecord getRecord(UUID sessionId)

}
