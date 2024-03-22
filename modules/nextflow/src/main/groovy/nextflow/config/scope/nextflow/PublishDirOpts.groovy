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

package nextflow.config.scope.nextflow

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString

/**
 * Model publishDir options
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@ToString(includePackage = false, includeNames = true)
@EqualsAndHashCode
@CompileStatic
class PublishDirOpts {

    static final public PublishDirOpts EMPTY = new PublishDirOpts(Map.of())


    final String mode
    final Boolean enabled
    final Boolean failOnError
    final String pattern
    final Object contentType
    final Boolean overwrite
    final String storageClass
    final Map<String,String> tags

    PublishDirOpts(Map opts) {
        mode = opts.mode
        enabled = asBool(opts.enabled)
        failOnError = asBool(opts.failOnError)
        overwrite = asBool(opts.overwrite)
        pattern = opts.pattern
        contentType = opts.contentType
        storageClass = opts.storageClass
        tags = opts.tags as Map<String,String>
    }

    private Boolean asBool(Object value) {
        if( value==null )
            return null
        return Boolean.valueOf(value as String)
    }
}
