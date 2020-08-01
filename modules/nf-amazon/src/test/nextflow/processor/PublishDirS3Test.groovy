/*
 * Copyright 2020, Seqera Labs
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

package nextflow.processor

import spock.lang.Specification

import java.nio.file.FileSystems

import nextflow.file.FileHelper

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class PublishDirS3Test extends Specification {

    def 'should change mode to `copy`' () {

        given:
        def processor = [:] as TaskProcessor
        processor.name = 'foo'

        def targetDir = FileHelper.asPath( 's3://bucket/work' )
        def publisher = new PublishDir(mode:'symlink', path: targetDir, sourceFileSystem: FileSystems.default)

        when:
        publisher.validatePublishMode()
        then:
        publisher.mode == PublishDir.Mode.COPY
    }


}
