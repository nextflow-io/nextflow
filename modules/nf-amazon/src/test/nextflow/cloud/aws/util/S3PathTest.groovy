package nextflow.cloud.aws.util

import nextflow.file.FileHelper
import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class S3PathTest extends Specification {

    @Unroll
    def 'should convert to uri string' () {

        expect:
        FileHelper.asPath(PATH).toUriString() == STR

        where:
        _ | PATH                | STR
        _ | 's3://foo'          | 's3://foo/'
        _ | 's3://foo/bar'      | 's3://foo/bar'
        _ | 's3://foo/b a r'    | 's3://foo/b a r'
        _ | 's3://f o o/bar'    | 's3://f o o/bar'
        _ | 's3://f_o_o/bar'    | 's3://f_o_o/bar'

    }

    @Unroll
    def 'should convert to string' () {

        expect:
        FileHelper.asPath(PATH).toString() == STR

        where:
        _ | PATH                | STR
        _ | 's3://foo'          | '/foo/'
        _ | 's3://foo/bar'      | '/foo/bar'
        _ | 's3://foo/b a r'    | '/foo/b a r'
        _ | 's3://f o o/bar'    | '/f o o/bar'
        _ | 's3://f_o_o/bar'    | '/f_o_o/bar'

    }

}
