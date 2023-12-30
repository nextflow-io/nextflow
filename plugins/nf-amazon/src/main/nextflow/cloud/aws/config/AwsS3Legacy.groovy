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

package nextflow.cloud.aws.config


import groovy.transform.CompileStatic
import nextflow.util.Duration
import nextflow.util.MemoryUnit
import org.apache.commons.lang.StringUtils
/**
 * Handle AWS S3 client legacy configuration
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class AwsS3Legacy {

    private Map config

    AwsS3Legacy( Map config ) {
        this.config = config
    }

    Map<String,?> getAwsClientConfig() {
        return config != null
                ? normalizeAwsClientConfig(config)
                : new HashMap<String,?>()
    }

    static protected Map normalizeAwsClientConfig(Map<String,?> client) {

        normalizeMemUnit(client, 'uploadChunkSize');
        normalizeDuration(client, 'uploadRetrySleep');


        def result = [:]
        client.each { String name, value ->
            def newKey = name.isCamelCase() ? StringUtils.splitByCharacterTypeCamelCase(name).join('_').toLowerCase() : name
            result.put(newKey,value?.toString())
        }
        return result
    }

    static void normalizeMemUnit(Map<String,?> client, String key) {
        if( client.get(key) instanceof String ) {
            client.put(key, MemoryUnit.of((String)client.get(key)))
        }
        if( client.get(key) instanceof MemoryUnit ) {
            client.put(key, ((MemoryUnit)client.get(key)).toBytes())
        }
    }

    static void normalizeDuration(Map<String,?> client, String key)  {
        if( client.get(key) instanceof String ) {
            client.put(key, Duration.of((String)client.get(key)))
        }
        if( client.get(key) instanceof Duration ) {
            client.put(key, ((Duration)client.get(key)).toMillis())
        }
    }
}
