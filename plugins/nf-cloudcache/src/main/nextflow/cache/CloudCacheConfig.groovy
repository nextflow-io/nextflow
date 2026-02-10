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

import groovy.transform.CompileStatic
import nextflow.config.spec.ConfigOption
import nextflow.config.spec.ConfigScope
import nextflow.config.spec.ScopeName
import nextflow.script.dsl.Description

/**
 * Model cloud cache configuration
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@ScopeName("cloudcache")
@Description("""
    The `cloudcache` scope controls the use of object storage as cache storage for workflow execution metadata.
""")
@CompileStatic
class CloudCacheConfig implements ConfigScope {

    @ConfigOption
    @Description("""
        Enable the use of cloud cache (default: `false`).
    """)
    final boolean enabled

    @ConfigOption
    @Description("""
        The path to the cloud storage bucket for cache metadata (e.g. `s3://bucket/cache`).
    """)
    final String path

    /* required by extension point -- do not remove */
    CloudCacheConfig() {}

    CloudCacheConfig(Map opts) {
        this.enabled = opts.enabled as boolean
        this.path = opts.path as String
    }

}
