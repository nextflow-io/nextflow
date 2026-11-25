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

package nextflow.file.http

import spock.lang.Specification

import nextflow.SysEnv

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class XFileSystemConfigTest extends Specification {

    def 'should create with default config settings' () {
        when:
        def config = new XFileSystemConfig()
        then:
        config.retryCodes() == XFileSystemConfig.DEFAULT_RETRY_CODES.tokenize(',').collect( it -> it as int )
        config.backOffDelay() == XFileSystemConfig.DEFAULT_BACK_OFF_DELAY
        config.backOffBase() == XFileSystemConfig.DEFAULT_BACK_OFF_BASE
        config.maxAttempts() == XFileSystemConfig.DEFAULT_MAX_ATTEMPTS
    }

    def 'should create with custom config settings' () {
        given:
        SysEnv.push([NXF_HTTPFS_MAX_ATTEMPTS: '10',
                     NXF_HTTPFS_BACKOFF_BASE: '300',
                     NXF_HTTPFS_DELAY       : '400',
                     NXF_HTTPFS_RETRY_CODES : '1,2,3'])

        when:
        def config = new XFileSystemConfig()
        then:
        config.retryCodes() == [1,2,3]
        config.backOffDelay() == 400
        config.backOffBase() == 300
        config.maxAttempts() == 10

        cleanup:
        SysEnv.pop()
    }

    def 'should create from an explicit environment' () {
        when:
        def config = new XFileSystemConfig([NXF_HTTPFS_MAX_ATTEMPTS: '5', NXF_HTTPFS_DELAY: '0'])
        then:
        config.maxAttempts() == 5
        config.backOffDelay() == 0
        config.backOffBase() == XFileSystemConfig.DEFAULT_BACK_OFF_BASE
        config.retryCodes() == XFileSystemConfig.DEFAULT_RETRY_CODES.tokenize(',').collect( it -> it as int )
    }

    def 'should replace and reload the global config' () {
        given:
        def custom = new XFileSystemConfig([NXF_HTTPFS_MAX_ATTEMPTS: '7'])

        when:
        XFileSystemConfig.reset(custom)
        then:
        XFileSystemConfig.config().is(custom)

        when:
        SysEnv.push([NXF_HTTPFS_MAX_ATTEMPTS: '9'])
        XFileSystemConfig.reset()
        then:
        XFileSystemConfig.config().maxAttempts() == 9

        cleanup:
        SysEnv.pop()
        XFileSystemConfig.reset()
    }
}
