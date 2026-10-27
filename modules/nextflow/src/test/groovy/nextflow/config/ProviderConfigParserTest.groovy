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

package nextflow.config

import java.nio.file.Files

import nextflow.exception.AbortOperationException
import nextflow.scm.ProviderConfig
import spock.lang.Specification

class ProviderConfigParserTest extends Specification {

    def 'should parse the scm config with environment variables' () {
        when:
        def config = ProviderConfigParser.parse('''
            providers {
                custom {
                    platform = 'gitlab'
                    token = env('HOME')
                }
            }
            ''')
        def providers = ProviderConfig.createFromMap(config)
        then:
        providers.find { it.name == 'custom' }.token == System.getenv('HOME')
        providers.find { it.name == 'custom' }.platform == 'gitlab'
    }

    def 'should return an empty config when the default scm file is missing' () {
        given:
        def folder = Files.createTempDirectory('test')
        def defaultFile = ProviderConfig.DEFAULT_SCM_FILE
        ProviderConfig.DEFAULT_SCM_FILE = folder.resolve('scm')

        expect:
        ProviderConfigParser.getFromFile(folder.resolve('scm')) == [:]

        cleanup:
        ProviderConfig.DEFAULT_SCM_FILE = defaultFile
        folder?.deleteDir()
    }

    def 'should fail when a custom scm file is missing' () {
        given:
        def folder = Files.createTempDirectory('test')

        when:
        ProviderConfigParser.getFromFile(folder.resolve('missing'))
        then:
        thrown(AbortOperationException)

        cleanup:
        folder?.deleteDir()
    }
}
