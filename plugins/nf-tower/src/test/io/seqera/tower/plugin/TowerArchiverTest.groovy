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

package io.seqera.tower.plugin

import java.nio.file.Files
import java.nio.file.Path
import java.util.concurrent.Executors

import nextflow.Session
import nextflow.file.FileHelper
import spock.lang.Specification
import spock.lang.Timeout
import spock.lang.Unroll
import test.TestHelper
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Timeout(5)
class TowerArchiverTest extends Specification {

    def 'should empty target path' () {
        expect:
        TowerArchiver.create(Mock(Session), [:]) == null
    }

    @Unroll
    def 'should final target path' () {
        given:
        def sess = Mock(Session) {
            getConfig() >> [:]
        }
        def helper = TowerArchiver.create(sess, [NXF_ARCHIVE_DIR:'/data,http://bucket/data/export'])

        expect:
        helper.getBaseDir() == Path.of('/data')
        helper.getTargetDir() == FileHelper.asPath('http://bucket/data/export')

        and:
        helper.archivePath(Path.of(SOURCE)) == (TARGET ? FileHelper.asPath(TARGET) : null)

        where:
        SOURCE                      | TARGET
        '/some/file.txt'            | null
        '/data/work/some/file.txt'  | 'http://bucket/data/export/work/some/file.txt'

    }

    def 'should parse paths' () {
        expect:
        TowerArchiver.parse(null) == []
        TowerArchiver.parse('/this,http://that') == ['/this','http://that']
        TowerArchiver.parse('/this,http://that/and/1\\,2\\,3') == ['/this','http://that/and/1,2,3']
    }

    def 'should failed parsing paths' () {
        when:
        TowerArchiver.parse('/this')
        then:
        thrown(IllegalArgumentException)

        when:
        TowerArchiver.parse('this,/that')
        then:
        thrown(IllegalArgumentException)

        when:
        TowerArchiver.parse('/this,that')
        then:
        thrown(IllegalArgumentException)

        when:
        TowerArchiver.parse('/a,/b,/c')
        then:
        thrown(IllegalArgumentException)
    }

    def 'should validate retry policy' () {
        given:
        def sess = Mock(Session) { getConfig() >> [:] }
        def archiver = new TowerArchiver(Path.of('/foo'),Path.of('/bar'), sess, [:])

        expect:
        !archiver.retryPattern().matcher('hello').find()
        !archiver.retryPattern().matcher('slow').find()
        and:
        archiver.retryPattern().matcher('slowdown').find()
        archiver.retryPattern().matcher('SlowDown').find()
        archiver.retryPattern().matcher('Slow Down').find()
        archiver.retryPattern().matcher('TooMany').find()
        archiver.retryPattern().matcher('Too Many').find()
    }

    def 'should archive path' () {
        given:
        def folder = TestHelper.createInMemTempDir()
        def source = folder.resolve('foo.txt'); source.text = 'xya'
        and:
        def sess = Mock(Session) { getConfig() >> [:] }
        def exec = Executors.newSingleThreadExecutor()
        def archiver = Spy(new TowerArchiver(folder,Path.of('/bar'), sess, [:]))
        def target = folder.resolve('bar.txt')

        when:
        archiver.archiveFile(source)
        then:
        1 * archiver.archivePath(source) >> target
        and:sleep 100   // <-- sleep a bit to allow the task to enter in the exec pool
        exec.shutdown()
        and:
        target.exists()
    }

    def 'should archive Nextflow log file' () {
        given:
        def workDir = Files.createTempDirectory('workdir')
        def targetDir = Files.createTempDirectory('target')
        def nextflowLog = workDir.resolve('nextflow.log'); nextflowLog.text = 'nextflow logs'
        def targetLog = targetDir.resolve('nextflow.log')
        and:
        def sess = Mock(Session) { getConfig() >> [:] }
        def exec = Executors.newSingleThreadExecutor()
        def archiver = new TowerArchiver(workDir,targetDir, sess, ['NXF_WORK': workDir.toString(), 'NXF_LOG_FILE': 'nextflow.log'])

        when:
        archiver.archiveLogs()
        sleep 100 // <-- sleep a bit to allow the task to enter in the exec pool
        exec.shutdown()

        then:
        targetLog.exists()

        cleanup:
        workDir?.deleteDir()
        targetDir?.deleteDir()
    }

}
