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

package io.seqera.tower.plugin

import java.nio.file.Files

import nextflow.SysEnv
import spock.lang.Specification
import test.TestHelper

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CacheCommandTest extends Specification {

    def cleanup() {
        SysEnv.pop()
    }

    def 'should backup logs to the work dir given by the config when cloud cache is enabled' () {
        given:
        def folder = Files.createTempDirectory('test')
        def remote = TestHelper.createInMemTempDir()
        def outFile = folder.resolve('nf-out.txt'); outFile.text = 'out file'
        def logFile = folder.resolve('nf-log.txt'); logFile.text = 'log file'
        and:
        SysEnv.push([NXF_CLOUDCACHE_PATH: 's3://bucket/cache', NXF_OUT_FILE: outFile.toString(), NXF_LOG_FILE: logFile.toString()])

        when:
        new CacheCommand().exec('nf-tower', 'cache-backup', [], [workDir: remote])
        then:
        remote.resolve('nf-out.txt').text == 'out file'
        remote.resolve('nf-log.txt').text == 'log file'

        cleanup:
        folder?.deleteDir()
    }

    def 'should backup logs to the work dir given by the env when cloud cache is not enabled' () {
        given:
        def folder = Files.createTempDirectory('test')
        def remote = folder.resolve('remote'); remote.mkdir()
        def outFile = folder.resolve('nf-out.txt'); outFile.text = 'out file'
        and:
        SysEnv.push([NXF_WORK: remote.toString(), NXF_UUID: UUID.randomUUID().toString(), NXF_OUT_FILE: outFile.toString()])

        when:
        new CacheCommand().exec('nf-tower', 'cache-backup', [], [workDir: folder.resolve('other').toString()])
        then:
        remote.resolve('nf-out.txt').text == 'out file'

        cleanup:
        folder?.deleteDir()
    }

}
