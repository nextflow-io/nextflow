/*
 * Copyright 2013-2024, Seqera Labs
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

import java.nio.file.Path

import com.google.common.hash.HashCode
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.SysEnv
import nextflow.cache.CacheDB
import nextflow.cache.CacheFactory
import nextflow.exception.AbortOperationException
import nextflow.fusion.FusionHelper
import nextflow.file.FileHelper
import nextflow.trace.TraceRecord
import nextflow.util.HistoryFile

import static nextflow.fusion.FusionConfig.FUSION_PATH

import static nextflow.util.StringUtils.*

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
        if( !args )
            throw new AbortOperationException("Missing id of the task to debug or remote task path")

        final criteria = args.pop()

        if ( isRemotePath(criteria) && args ) {
            final image = args.pop()
            runRemoteTask(criteria, image)
        }
        else {
            runCacheTask(criteria)
        }

    }

    protected runRemoteTask(String workDir, String image) {
        final path = FileHelper.asPath(workDir)
        final cmd = buildWaveRunCmd(path.scheme)
        if ( image.startsWith("wave.seqera.io/") ) {
            // Run a waved container
            cmd.runContainer(image, buildCommand(path))
        }
        else {
            // Wave it before running
            List<String> args = [image] + buildCommand(path)
            cmd.runContainer(args)
        }
    }

    protected runCacheTask(String criteria) {
        HistoryFile history = HistoryFile.DEFAULT
        HistoryFile.Record runRecord
        if( !history.exists() || history.empty() || !(runRecord=history.findByIdOrName('last')[0]) )
            throw new AbortOperationException("It looks no pipeline was executed in this folder (or execution history is empty)")

        this.cacheDb = CacheFactory.create(runRecord.sessionId, runRecord.runName)
        cacheDb.openForRead()
        try {
            final trace = getOrFindTrace(criteria)
            if( trace==null )
                throw new AbortOperationException("Cannot find any task with id: '$criteria'")
            if( !isRemotePath(trace.workDir) )
                throw new AbortOperationException("Cannot run non-fusion enabled task - Task work dir: $trace.workDir")
            log.info "Launching debug session for task '${trace.get('name')}' - work directory: ${trace.workDir}"
            final cmd = buildWaveRunCmd(getUrlProtocol(trace.workDir))
            cmd.runContainer(trace.get('container')?.toString(), buildCommand(trace.workDir))
        }
        finally {
            cacheDb.close()
        }
    }

    protected static List<String> buildCommand(String workDir) {
        final workPath = FileHelper.asPath(workDir)
        return buildCommand(workPath)
    }

    protected static List<String> buildCommand(Path workPath) {
        final fusionPath = FusionHelper.toContainerMount(workPath, workPath.scheme)
        return [FUSION_PATH, 'sh', '-c', "\'cd $fusionPath && PS1=\"[fusion] \" exec bash --norc\'".toString()]
    }

    protected WaveRunCmd buildWaveRunCmd(String scheme) {
        final result = new WaveRunCmd(session)
        result.withContainerParams([tty:true, privileged: true])
        if( scheme=='s3' ) {
            result.withEnvironment('AWS_ACCESS_KEY_ID')
            result.withEnvironment('AWS_SECRET_ACCESS_KEY')
        }
        else if( scheme=='gs' && SysEnv.containsKey('GOOGLE_APPLICATION_CREDENTIALS') ) {
            final path = Path.of(SysEnv.get('GOOGLE_APPLICATION_CREDENTIALS'))
            result.withMounts(List.of(path))
            result.withEnvironment('GOOGLE_APPLICATION_CREDENTIALS')
        }
        return result
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

    static boolean isRemotePath(String path) {
        if (!path) return false
        final result = getUrlProtocol(path)
        return result!=null && result!='file'
    }
}
