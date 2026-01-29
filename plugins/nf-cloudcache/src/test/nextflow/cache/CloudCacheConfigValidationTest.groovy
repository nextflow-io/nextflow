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

package nextflow.cache

import nextflow.config.ConfigValidator
import nextflow.plugin.Plugins
import org.junit.Rule
import spock.lang.Specification
import test.OutputCapture

/**
 * Integration test to verify cloudcache config options are recognized by ConfigValidator.
 * This replicates the issue reported in https://github.com/nextflow-io/nextflow/issues/6773
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CloudCacheConfigValidationTest extends Specification {

    @Rule
    OutputCapture capture = new OutputCapture()

    def setupSpec() {
        // Start the plugin system and load the cloudcache plugin
        Plugins.init()
        Plugins.start('nf-cloudcache')
    }

    def cleanupSpec() {
        Plugins.stop()
    }

    def 'should not warn about cloudcache.enabled config option' () {
        when:
        new ConfigValidator().validate([
            cloudcache: [
                enabled: true
            ]
        ])
        then:
        !capture.toString().contains("Unrecognized config option 'cloudcache.enabled'")
    }

    def 'should not warn about cloudcache.path config option' () {
        when:
        new ConfigValidator().validate([
            cloudcache: [
                path: 's3://bucket/cache'
            ]
        ])
        then:
        !capture.toString().contains("Unrecognized config option 'cloudcache.path'")
    }

    def 'should not warn about any cloudcache config options' () {
        when:
        new ConfigValidator().validate([
            cloudcache: [
                enabled: true,
                path: 's3://bucket/cache'
            ]
        ])
        then:
        !capture.toString().contains("Unrecognized config option 'cloudcache")
    }

    def 'should still warn about invalid cloudcache config options' () {
        when:
        new ConfigValidator().validate([
            cloudcache: [
                enabled: true,
                path: 's3://bucket/cache',
                invalidOption: 'foo'
            ]
        ])
        then:
        capture.toString().contains("Unrecognized config option 'cloudcache.invalidOption'")
        !capture.toString().contains("Unrecognized config option 'cloudcache.enabled'")
        !capture.toString().contains("Unrecognized config option 'cloudcache.path'")
    }

}
