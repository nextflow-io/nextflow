package nextflow.cloud.aws.util

import nextflow.cloud.aws.nio.S3Path
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

    def 'should check equals and hashcode' () {
        given:
        def path1 = FileHelper.asPath('s3://foo/some/foo.txt')
        def path2 = FileHelper.asPath('s3://foo/some/foo.txt')
        def path3 = FileHelper.asPath('s3://foo/some/bar.txt')
        def path4 = FileHelper.asPath('s3://bar/some/foo.txt')

        expect:
        path1 == path2
        path1 != path3
        path3 != path4
        and:
        path1.hashCode() == path2.hashCode()
        path1.hashCode() != path3.hashCode()
        path3.hashCode() != path4.hashCode()
    }

    @Unroll
    def 'should determine bucket name' () {
        expect:
        S3Path.bucketName(new URI(URI_PATH)) == BUCKET

        where:
        URI_PATH            | BUCKET
        's3:///'            | null
        's3:///foo'         | 'foo'
        's3:///foo/'        | 'foo'
        's3:///foo/bar'     | 'foo'
    }

    @Unroll
    def 'should normalise path' () {
        expect:
        FileHelper.asPath(PATH).normalize() == FileHelper.asPath(EXPECTED)

        where:
        PATH                        | EXPECTED
        's3://foo'                  | 's3://foo'
        's3://foo/x/y/z.txt'        | 's3://foo/x/y/z.txt'
        's3://foo/x/y/./z.txt'      | 's3://foo/x/y/z.txt'
        's3://foo/x/y/../z.txt'     | 's3://foo/x/z.txt'
        's3://foo/x/y/../../z.txt'  | 's3://foo/z.txt'
        's3://foo/x/y//z.txt'       | 's3://foo/x/y/z.txt'
        's3://foo/./z.txt'          | 's3://foo/z.txt'
    }
}
