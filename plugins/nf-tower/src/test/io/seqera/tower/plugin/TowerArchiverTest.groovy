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

import java.nio.file.Path

import nextflow.Session
import nextflow.file.FileHelper
import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
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

}
