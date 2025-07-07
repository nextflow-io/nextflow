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

package io.seqera.tower.plugin

import static java.nio.file.StandardCopyOption.*

import java.nio.file.FileSystems
import java.nio.file.Path
import java.nio.file.Paths

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.exception.AbortOperationException
import nextflow.file.FileHelper
/**
 * Backup Nextflow logs, timeline and reports files
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class LogsHandler {

    @PackageScope Path localOutFile
    @PackageScope Path localLogFile
    @PackageScope Path localTimelineFile
    @PackageScope Path localTowerConfig
    @PackageScope Path localTowerReports
    @PackageScope Path remoteWorkDir

    @PackageScope Path getRemoteOutFile() { remoteWorkDir.resolve(localOutFile.getName()) }
    @PackageScope Path getRemoteLogFile() { remoteWorkDir.resolve(localLogFile.getName()) }
    @PackageScope Path getRemoteTimelineFile() { remoteWorkDir.resolve(localTimelineFile.getName()) }
    @PackageScope Path getRemoteTowerConfig() { remoteWorkDir.resolve(localTowerConfig.getName()) }
    @PackageScope Path getRemoteTowerReports() { remoteWorkDir.resolve(localTowerReports.getName()) }

    LogsHandler(Session session, Map<String,String> env) {
        if( !session.workDir )
            throw new AbortOperationException("Missing workflow work directory")
        if( session.workDir.fileSystem == FileSystems.default )
            throw new AbortOperationException("Logs handler is only meant to be used with a remote workflow work directory")
        this.remoteWorkDir = session.workDir

        if( env.NXF_OUT_FILE )
            localOutFile = Paths.get(env.NXF_OUT_FILE)
        if( env.NXF_LOG_FILE )
            localLogFile = Paths.get(env.NXF_LOG_FILE)
        if( env.NXF_TML_FILE )
            localTimelineFile = Paths.get(env.NXF_TML_FILE)
        if( env.TOWER_CONFIG_FILE )
            localTowerConfig = Paths.get(env.TOWER_CONFIG_FILE)
        if( env.TOWER_REPORTS_FILE )
            localTowerReports = Paths.get(env.TOWER_REPORTS_FILE)
    }

    void saveFiles() {
        log.trace "Checkpointing logs, timeline and report files"
        // — upload out file
        try {
            if( localOutFile?.exists() )
                FileHelper.copyPath(localOutFile, remoteOutFile, REPLACE_EXISTING)
        }
        catch (Throwable e) {
            log.warn "Unable to upload nextflow out file: $localOutFile — reason: ${e.message ?: e}", e
        }
        // — upload log file
        try {
            if( localLogFile?.exists() )
                FileHelper.copyPath(localLogFile, remoteLogFile, REPLACE_EXISTING)
        }
        catch (Throwable e) {
            log.warn "Unable to upload nextflow log file: $localLogFile — reason: ${e.message ?: e}", e
        }
        // — upload timeline file
        try {
            if( localTimelineFile?.exists() )
                FileHelper.copyPath(localTimelineFile, remoteTimelineFile, REPLACE_EXISTING)
        }
        catch (Throwable e) {
            log.warn "Unable to upload nextflow timeline file: $localTimelineFile — reason: ${e.message ?: e}", e
        }
        // — upload tower config file
        try {
            if( localTowerConfig?.exists() )
                FileHelper.copyPath(localTowerConfig, remoteTowerConfig, REPLACE_EXISTING)
        }
        catch (Throwable e) {
            log.warn "Unable to upload tower config file: $localTowerConfig — reason: ${e.message ?: e}", e
        }
        // — upload tower reports file
        try {
            if( localTowerReports?.exists() )
                FileHelper.copyPath(localTowerReports, remoteTowerReports, REPLACE_EXISTING)
        }
        catch (Throwable e) {
            log.warn "Unable to upload tower reprts file: $localTowerReports — reason: ${e.message ?: e}", e
        }
    }

}
