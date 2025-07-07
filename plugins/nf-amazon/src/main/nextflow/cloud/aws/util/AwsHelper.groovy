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

package nextflow.cloud.aws.util

import software.amazon.awssdk.services.s3.model.ObjectCannedACL
import com.google.common.base.CaseFormat

/**
 * Helper class for AWS
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AwsHelper {

    static ObjectCannedACL parseS3Acl(String value) {
        if( !value )
            return null

        return value.contains('-')
                ? ObjectCannedACL.valueOf(CaseFormat.LOWER_HYPHEN.to(CaseFormat.UPPER_UNDERSCORE, value))
                : ObjectCannedACL.valueOf(CaseFormat.UPPER_CAMEL.to(CaseFormat.UPPER_UNDERSCORE,value))
    }

}
