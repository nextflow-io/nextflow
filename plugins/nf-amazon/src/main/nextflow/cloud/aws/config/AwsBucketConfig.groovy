/*
 * Copyright 2020-2025, Seqera Labs
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
import groovy.util.logging.Slf4j
import nextflow.config.spec.ConfigOption
import nextflow.script.dsl.Description

/**
 * Model AWS S3 bucket config settings
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
@CompileStatic
class AwsBucketConfig extends AwsS3CommonConfig {

    @ConfigOption
    @Description("""
        S3 Bucket specific AWS region (e.g. `us-east-1`).
    """)
    final String region

    AwsBucketConfig(Map opts) {
        super(opts)
        this.region = opts.region as String
    }

    Map<String, Object> toBucketConfigMap(){
        final map = super.toBucketConfigMap()
        if( region ) map.put('region', region)
        return map
    }
}
