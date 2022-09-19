/*
 * Copyright 2020-2022, Seqera Labs
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
 *
 */

package io.seqera.wave.plugin.cli

import com.google.common.hash.HashCode
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.cache.CacheDB
import nextflow.cache.CacheFactory
import nextflow.exception.AbortOperationException
import nextflow.executor.fusion.FusionHelper
import nextflow.file.FileHelper
import nextflow.trace.TraceRecord
import nextflow.util.HistoryFile

/**
 * Implement wave container task debug
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class WaveDebugCmd {

    private CacheDB cacheDb

    private Session session

    WaveDebugCmd(Session session) {
        this.session = session
    }

    protected void taskDebug(List<String> args) {
        HistoryFile history = HistoryFile.DEFAULT
        HistoryFile.Record runRecord
        if( !history.exists() || history.empty() || !(runRecord=history.findByIdOrName('last')[0]) )
            throw new AbortOperationException("It looks no pipeline was executed in this folder (or execution history is empty)")

        if( !args )
            throw new AbortOperationException("Missing id of the task to debug")
        final criteria = args.pop()

        this.cacheDb = CacheFactory.create(runRecord.sessionId, runRecord.runName)
        cacheDb.openForRead()
        try {
            final trace = getOrFindTrace(criteria)
            if( trace==null )
                throw new AbortOperationException("Cannot find any task with id: '$criteria'")

            if( !trace.workDir.startsWith('s3://') )
                throw new AbortOperationException("Cannot run non-fusion enabled task - Task work dir: $trace.workDir")

            log.info "Launching debug session for task '${trace.get('name')}' - work directory: ${trace.workDir}"
            final workDir = FileHelper.asPath(trace.workDir)
            final fusionPath = FusionHelper.toContainerMount(workDir, 's3')
            new WaveRunCmd(session)
                    .withContainerParams([tty:true, privileged: true])
                    .withEnvironment('AWS_ACCESS_KEY_ID')
                    .withEnvironment('AWS_SECRET_ACCESS_KEY')
                    .withEnvironment("NXF_FUSION_WORK=$fusionPath".toString())
                    .runContainer(trace.get('container')?.toString(), ['sh', '-c', '\'cd $NXF_FUSION_WORK && exec bash\''])
        }
        finally {
            cacheDb.close()
        }
    }

    protected TraceRecord getOrFindTrace(String criteria) {
        TraceRecord result = null
        final norm = criteria.replace('/','')
        if( norm.size()==32 ) try {
            final hash = HashCode.fromString(norm)
            result = cacheDb.getTraceRecord(hash)
        }
        catch (IllegalArgumentException e) {
            log.debug "Not a valid task hash: $criteria"
        }

        if( !result ) {
            final condition = { TraceRecord trace ->
                final hash = trace.get('hash') as String
                final name = trace.get('name') as String
                return hash?.contains(criteria)
                        || hash.replace('/','')?.contains(criteria)
                        || name?.contains(criteria)
                        || trace.workDir?.contains(criteria)
            }
            result = cacheDb.findTraceRecord(condition)
        }

        return result
    }

}
