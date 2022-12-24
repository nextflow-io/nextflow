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

package nextflow.cloud.aws.fusion

import groovy.transform.CompileStatic
import nextflow.Global
import nextflow.fusion.FusionConfig
import nextflow.fusion.FusionEnv
import org.pf4j.Extension
/**
 * Implements {@link FusionEnv} for AWS cloud
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Extension
@CompileStatic
class AwsFusionEnv implements FusionEnv {

    @Override
    Map<String, String> getEnvironment(String scheme, FusionConfig config) {
        if( scheme!='s3' )
            return Collections.<String,String>emptyMap()

        final result = new HashMap<String,String>()
        final endpoint = Global.getAwsS3Endpoint()
        final creds = config.exportAwsAccessKeys() ? Global.getAwsCredentials() : Collections.<String>emptyList()
        if( creds ) {
            result.AWS_ACCESS_KEY_ID = creds[0]
            result.AWS_SECRET_ACCESS_KEY = creds[1]
        }
        if( endpoint )
            result.AWS_S3_ENDPOINT = endpoint
        return result
    }
}
