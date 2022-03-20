/*
 * Copyright (c) 2019, Seqera Labs.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * This Source Code Form is "Incompatible With Secondary Licenses", as
 * defined by the Mozilla Public License, v. 2.0.
 */

package io.seqera.tower.plugin

import static java.nio.file.StandardCopyOption.*

import java.nio.file.NoSuchFileException
import java.nio.file.Path
import java.nio.file.Paths

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.extension.FilesEx
import nextflow.file.FileHelper
import nextflow.plugin.BasePlugin
import org.pf4j.PluginWrapper
/**
 * Nextflow Tower plugin
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class TowerPlugin extends BasePlugin {

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

    TowerPlugin(PluginWrapper wrapper) {
        super(wrapper)
        init()
    }


    protected void init(Map<String,String> env = System.getenv()) {
        if( !env.NXF_UUID )
            return
        if( !env.NXF_REMOTE_WORK )
            return

        sessionUuid = env.NXF_UUID
        localCachePath = Paths.get(".nextflow/cache/${sessionUuid}")
        remoteWorkDir = FileHelper.asPath(env.NXF_REMOTE_WORK)

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

    @Override
    void start() {
        super.start()
        restoreCacheFiles()
    }

    @Override
    void stop() {
        saveCacheFiles()
        super.stop()
    }

    protected void restoreCacheFiles() {
        if( !remoteWorkDir || !sessionUuid )
            return

        try {
            log.debug "Restoring Nextflow cache: $remoteCachePath ==> $localCachePath"
            localCachePath.deleteDir()
            localCachePath.parent.mkdirs()
            FileHelper.copyPath(remoteCachePath, localCachePath, REPLACE_EXISTING)
        }
        catch (NoSuchFileException e) {
            log.debug "Remote cache restore ignored -- reason: ${e.message ?: e}"
        }
    }

    protected void saveCacheFiles() {
        if( !remoteWorkDir || !sessionUuid )
            return
        
        // upload nextflow cache metadata
        try {
            log.debug "Saving Nextflow cache: $localCachePath ==> $remoteCachePath"
            remoteCachePath.deleteDir()
            remoteCachePath.parent.mkdirs()
            FilesEx.copyTo(localCachePath, remoteCachePath)
        }
        catch (Throwable e) {
            log.warn "Failed to backup resume metadata to remote store path: ${remoteCachePath.toUriString()}", e
        }
        // -- upload out file
        try {
            if( localOutFile?.exists() )
                FileHelper.copyPath(localOutFile, remoteOutFile, REPLACE_EXISTING)
        }
        catch (Throwable e) {
            log.warn "Unable to upload nextflow out file: $localOutFile -- reason: ${e.message ?: e}", e
        }
        // -- upload log file
        try {
            if( localLogFile?.exists() )
                FileHelper.copyPath(localLogFile, remoteLogFile, REPLACE_EXISTING)
        }
        catch (Throwable e) {
            log.warn "Unable to upload nextflow log file: $localLogFile -- reason: ${e.message ?: e}", e
        }
        // -- upload timeline file
        try {
            if( localTimelineFile?.exists() )
                FileHelper.copyPath(localTimelineFile, remoteTimelineFile, REPLACE_EXISTING)
        }
        catch (Throwable e) {
            log.warn "Unable to upload nextflow timeline file: $localTimelineFile -- reason: ${e.message ?: e}", e
        }
        // -- upload tower config file
        try {
            if( localTowerConfig?.exists() )
                FileHelper.copyPath(localTowerConfig, remoteTowerConfig, REPLACE_EXISTING)
        }
        catch (Throwable e) {
            log.warn "Unable to upload tower config file: $localTowerConfig -- reason: ${e.message ?: e}", e
        }
        // -- upload tower reports file
        try {
            if( localTowerReports?.exists() )
                FileHelper.copyPath(localTowerReports, remoteTowerReports, REPLACE_EXISTING)
        }
        catch (Throwable e) {
            log.warn "Unable to upload tower reprts file: $localTowerReports -- reason: ${e.message ?: e}", e
        }
    }
}
