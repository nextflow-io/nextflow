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
 */

package nextflow.cloud.aws.fusion

import groovy.transform.CompileStatic
import nextflow.SysEnv
import nextflow.cloud.aws.config.AwsConfig
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
        final awsConfig = AwsConfig.config()
        final endpoint = awsConfig.s3Config.endpoint
        final creds = config.exportStorageCredentials() ? awsCreds(awsConfig) : List.<String>of()
        if( creds ) {
            result.AWS_ACCESS_KEY_ID = creds[0]
            result.AWS_SECRET_ACCESS_KEY = creds[1]

            if( creds.size() > 2 )
                result.AWS_SESSION_TOKEN = creds[2]
        }
        if( endpoint )
            result.AWS_S3_ENDPOINT = endpoint
        if( awsConfig.region && awsConfig.s3Config.isCustomEndpoint() )
            result.FUSION_AWS_REGION = awsConfig.region
        if( awsConfig.s3Config.storageEncryption )
            result.FUSION_AWS_SERVER_SIDE_ENCRYPTION = awsConfig.s3Config.storageEncryption
        if( awsConfig.s3Config.storageKmsKeyId )
            result.FUSION_AWS_SSEKMS_KEY_ID = awsConfig.s3Config.storageKmsKeyId
        return result
    }

    protected List<String> awsCreds(AwsConfig awsConfig) {
        final result = awsConfig.getCredentials()
        if( result )
            return result

        if( SysEnv.get('AWS_ACCESS_KEY_ID') && SysEnv.get('AWS_SECRET_ACCESS_KEY') && SysEnv.get('AWS_SESSION_TOKEN') )
            return List.<String>of(SysEnv.get('AWS_ACCESS_KEY_ID'), SysEnv.get('AWS_SECRET_ACCESS_KEY'), SysEnv.get('AWS_SESSION_TOKEN'))

        if( SysEnv.get('AWS_ACCESS_KEY_ID') && SysEnv.get('AWS_SECRET_ACCESS_KEY') )
            return List.<String>of(SysEnv.get('AWS_ACCESS_KEY_ID'), SysEnv.get('AWS_SECRET_ACCESS_KEY'))
        else
            return List.<String>of()
    }
}
