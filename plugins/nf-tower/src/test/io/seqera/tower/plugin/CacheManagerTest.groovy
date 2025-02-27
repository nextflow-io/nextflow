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

import java.nio.file.Files
import java.nio.file.Paths

import nextflow.exception.AbortOperationException
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CacheManagerTest extends Specification {

    def 'should init empty files' () {
        when:
        new CacheManager([:])
        then:
        thrown(AbortOperationException)
    }

    def 'should upload cache files' () {
        given:
        def folder = Files.createTempDirectory('test')
        def remote = folder.resolve('remote'); remote.mkdir()
        def local = folder.resolve('local'); local.mkdir()
        def outFile = local.resolve('nf-out.txt'); outFile.text = 'out file'
        def logFile = local.resolve('nf-log.txt'); logFile.text = 'log file'
        def tmlFile = local.resolve('nf-tml.txt'); tmlFile.text = 'tml file'
        def cfgFile = local.resolve('tw-config.txt'); cfgFile.text = 'config file'
        def repFile = local.resolve('tw-report.txt'); repFile.text = 'report file'
        and:
        def uuid = UUID.randomUUID().toString()
        and:
        def ENV = [
                NXF_UUID:uuid,
                NXF_WORK: remote.toString(),
                NXF_OUT_FILE: outFile.toString(),
                NXF_LOG_FILE: logFile.toString(),
                NXF_TML_FILE: tmlFile.toString(),
                TOWER_CONFIG_FILE: cfgFile.toString(),
                TOWER_REPORTS_FILE: repFile.toString(),
        ]

        when:
        def tower = new CacheManager(ENV)
        then:
        tower.sessionUuid == uuid
        tower.localCachePath == Paths.get(".nextflow/cache/$uuid")
        tower.localOutFile == outFile
        tower.localLogFile == logFile
        tower.localTimelineFile == tmlFile
        tower.localTowerConfig == cfgFile
        tower.localTowerReports == repFile
        and:
        tower.remoteWorkDir == remote
        and:
        tower.remoteCachePath == remote.resolve(".nextflow/cache/$uuid")
        tower.remoteOutFile == remote.resolve( outFile.name )
        tower.remoteLogFile == remote.resolve( logFile.name )
        tower.remoteTimelineFile == remote.resolve( tmlFile.name )
        tower.remoteTowerConfig == remote.resolve( cfgFile.name )
        tower.remoteTowerReports == remote.resolve( repFile.name )

        when:
        // create local cache fake data
        tower.localCachePath = local.resolve(".nextflow/cache/$uuid");
        tower.localCachePath.mkdirs()
        tower.localCachePath.resolve('index-foo').text = 'index foo'
        tower.localCachePath.resolve('db').mkdir()
        tower.localCachePath.resolve('db/xxx').text = 'data xxx'
        tower.localCachePath.resolve('db/yyy').text = 'data yyy'
        and:
        tower.saveCacheFiles()
        then:
        tower.remoteCachePath.resolve('index-foo').text == 'index foo'
        tower.remoteCachePath.resolve('db/xxx').text == 'data xxx'
        tower.remoteCachePath.resolve('db/yyy').text == 'data yyy'
        and:
        tower.remoteOutFile.text == outFile.text
        tower.remoteLogFile.text == logFile.text
        tower.remoteTimelineFile.text == tmlFile.text
        tower.remoteTowerConfig.text == cfgFile.text
        tower.remoteTowerReports.text == repFile.text

        // simulate a 2nd run with different data
        when:
        tower.localCachePath.deleteDir()
        tower.localCachePath.mkdirs()
        tower.localCachePath.resolve('index-bar').text = 'index bar'
        tower.localCachePath.resolve('db').mkdir()
        tower.localCachePath.resolve('db/alpha').text = 'data alpha'
        tower.localCachePath.resolve('db/delta').text = 'data delta'
        and:
        tower.saveCacheFiles()
        then:
        tower.remoteCachePath.resolve('index-bar').text == 'index bar'
        tower.remoteCachePath.resolve('db/alpha').text == 'data alpha'
        tower.remoteCachePath.resolve('db/delta').text == 'data delta'
        and:
        !tower.remoteCachePath.resolve('index-foo').exists()
        !tower.remoteCachePath.resolve('db/xxx').exists()
        !tower.remoteCachePath.resolve('db/yyy').exists()
        and:
        tower.remoteOutFile.text == outFile.text
        tower.remoteLogFile.text == logFile.text
        tower.remoteTimelineFile.text == tmlFile.text
        tower.remoteTowerConfig.text == cfgFile.text
        tower.remoteTowerReports.text == repFile.text

        cleanup:
        folder?.deleteDir()
    }

    def 'should download cache files' () {
        given:
        def uuid = UUID.randomUUID().toString()
        def folder = Files.createTempDirectory('test')
        def local = folder.resolve('local'); local.mkdir()
        def outFile = local.resolve('nf-out.txt');
        def logFile = local.resolve('nf-log.txt')
        def tmlFile = local.resolve('nf-tml.txt')
        def cfgFile = local.resolve('tw-config.txt')
        def repFile = local.resolve('tw-report.txt')
        and:
        def remote = folder.resolve('remote'); remote.mkdir()
        remote.resolve('nf-out.txt').text = 'the out file'
        remote.resolve('nf-log.txt').text = 'the log file'
        remote.resolve('nf-tml.txt').text = 'the timeline file'
        remote.resolve('nf-config.txt').text = 'the config file'
        remote.resolve('nf-report.txt').text = 'the report file'
        and:
        remote.resolve(".nextflow/cache/$uuid").mkdirs()
        remote.resolve(".nextflow/cache/$uuid").resolve('index-bar').text = 'index bar'
        remote.resolve(".nextflow/cache/$uuid").resolve('db').mkdirs()
        remote.resolve(".nextflow/cache/$uuid").resolve('db/alpha').text = 'data alpha'
        remote.resolve(".nextflow/cache/$uuid").resolve('db/delta').text = 'data delta'
        and:
        def tower = new CacheManager([NXF_UUID: uuid, NXF_WORK: remote.toString()])

        when:
        tower.restoreCacheFiles()
        then:
        tower.localCachePath.resolve('index-bar').text == 'index bar'
        tower.localCachePath.resolve('db/alpha').text == 'data alpha'
        tower.localCachePath.resolve('db/delta').text == 'data delta'

        cleanup:
        folder?.deleteDir()
    }
}
