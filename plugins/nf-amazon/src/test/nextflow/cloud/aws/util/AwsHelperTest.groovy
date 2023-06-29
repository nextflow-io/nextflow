/*
 * Copyright (c) 2020-2021. Seqera Labs, S.L.
 *
 * All Rights reserved
 *
 */

package nextflow.cloud.aws.util

import com.amazonaws.services.s3.model.CannedAccessControlList
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AwsHelperTest extends Specification {

    def 'should parse S3 acl' () {
        expect:
        AwsHelper.parseS3Acl('PublicRead') == CannedAccessControlList.PublicRead
        AwsHelper.parseS3Acl('public-read') == CannedAccessControlList.PublicRead

        when:
        AwsHelper.parseS3Acl('unknown')
        then:
        thrown(IllegalArgumentException)
    }
}
