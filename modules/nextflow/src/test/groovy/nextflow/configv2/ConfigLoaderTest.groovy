/*
 * Copyright 2013-2023, Seqera Labs
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

package nextflow.configv2

import java.nio.file.Files
import java.nio.file.Path

import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ConfigLoaderTest extends Specification {

    def 'should load config file' () {
        given:
        def loader = new ConfigLoader()
        def folder = Files.createTempDirectory('test')
        and:
        def config = folder.resolve('nextflow.config.yml')
        config.text = '''
nextflow:
    workDir: '/some/path'
'''
        when:
        def nextflow = loader.load([config]).nextflowConfig()
        then:
        nextflow.workDir == Path.of('/some/path')

        cleanup:
        folder?.deleteDir()
    }
}
