/*
 * Copyright (c) 2020-2021. Seqera Labs, S.L.
 *
 * All Rights reserved
 *
 */

package nextflow.cloud.aws.util

import com.amazonaws.services.s3.model.CannedAccessControlList
import com.google.common.base.CaseFormat

/**
 * Helper class for AWS
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AwsHelper {

    static CannedAccessControlList parseS3Acl(String value) {
        if( !value )
            return null

        return value.contains('-')
                ? CannedAccessControlList.valueOf(CaseFormat.LOWER_HYPHEN.to(CaseFormat.UPPER_CAMEL,value))
                : CannedAccessControlList.valueOf(value)
    }

}
