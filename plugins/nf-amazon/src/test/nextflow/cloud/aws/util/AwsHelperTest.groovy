/*
 * Copyright (c) 2020-2021. Seqera Labs, S.L.
 *
 * All Rights reserved
 *
 */

package nextflow.cloud.aws.util

import software.amazon.awssdk.services.s3.model.ObjectCannedACL
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AwsHelperTest extends Specification {

    def 'should parse S3 acl' () {
        expect:
        AwsHelper.parseS3Acl('PublicRead') == ObjectCannedACL.PUBLIC_READ
        AwsHelper.parseS3Acl('public-read') == ObjectCannedACL.PUBLIC_READ
        AwsHelper.parseS3Acl('Private') == ObjectCannedACL.PRIVATE
        AwsHelper.parseS3Acl('private') == ObjectCannedACL.PRIVATE
        when:
        AwsHelper.parseS3Acl('unknown')
        then:
        thrown(IllegalArgumentException)
    }
}
