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

import java.nio.file.Files
import java.nio.file.NoSuchFileException
import java.nio.file.Path
import java.nio.file.Paths

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.Const
import nextflow.exception.AbortOperationException
import nextflow.extension.FilesEx
import nextflow.file.FileHelper
/**
 * Back and restore Nextflow cache
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class CacheManager {

    @PackageScope String sessionUuid
    @PackageScope Path localCachePath
    @PackageScope Path localOutFile
    @PackageScope Path localLogFile
    @PackageScope Path localTimelineFile
    @PackageScope Path localTowerConfig
    @PackageScope Path localTowerReports
    @PackageScope Path remoteWorkDir

    @PackageScope Path getRemoteCachePath() { remoteWorkDir.resolve(".nextflow/cache/${sessionUuid}") }
    @PackageScope Path getRemoteOutFile() { remoteWorkDir.resolve(localOutFile.getName()) }
    @PackageScope Path getRemoteLogFile() { remoteWorkDir.resolve(localLogFile.getName()) }
    @PackageScope Path getRemoteTimelineFile() { remoteWorkDir.resolve(localTimelineFile.getName()) }
    @PackageScope Path getRemoteTowerConfig() { remoteWorkDir.resolve(localTowerConfig.getName()) }
    @PackageScope Path getRemoteTowerReports() { remoteWorkDir.resolve(localTowerReports.getName()) }

    CacheManager(Map<String,String> env) {
        final work = env.get('NXF_WORK') ?: env.get('NXF_TEST_WORK')
        if( !work )
            throw new AbortOperationException("Missing target work dir - cache sync cannot be performed")
        this.remoteWorkDir = FileHelper.asPath(work)

        this.sessionUuid = env.get('NXF_UUID')
        if( !sessionUuid )
            throw new AbortOperationException("Missing target uuid - cache sync cannot be performed")

        this.localCachePath = Const.appCacheDir.resolve("cache/${sessionUuid}")

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

    protected void restoreCacheFiles() {
        if( !remoteWorkDir || !sessionUuid )
            return

        if(!Files.exists(remoteCachePath)) {
            log.debug "Remote cache path does not exist: $remoteCachePath - skipping cache restore"
            return
        }

        try {
            log.info "Restoring cache: ${remoteCachePath.toUriString()} => ${localCachePath.toUriString()}"
            localCachePath.deleteDir()
            localCachePath.parent.mkdirs()
            FileHelper.copyPath(remoteCachePath, localCachePath, REPLACE_EXISTING)
        }
        catch (NoSuchFileException e) {
            log.info "Remote cache restore ignored — reason: ${e.message ?: e}"
        }
    }

    protected void saveCacheFiles() {
        if( !remoteWorkDir || !sessionUuid )
            return

        if( !Files.exists(localCachePath) ) {
            log.debug "Local cache path does not exist: $localCachePath — skipping cache backup"
            return
        }

        // upload nextflow cache metadata
        try {
            log.info "Saving cache: ${localCachePath.toUriString()} => ${remoteCachePath.toUriString()}"
            remoteCachePath.deleteDir()
            remoteCachePath.parent.mkdirs()
            FilesEx.copyTo(localCachePath, remoteCachePath)
        }
        catch (Throwable e) {
            log.warn "Failed to backup resume metadata to remote store path: ${remoteCachePath.toUriString()} — cause: ${e}", e
        }

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
