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

package nextflow.cache

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CloudCacheConfigTest extends Specification {

    def 'should create empty config' () {
        when:
        def config = new CloudCacheConfig([:])
        then:
        !config.enabled
        config.path == null
    }

    def 'should create config with all options' () {
        when:
        def config = new CloudCacheConfig([enabled: true, path: 's3://bucket/cache'])
        then:
        config.enabled
        config.path == 's3://bucket/cache'
    }

}
