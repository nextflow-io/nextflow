package nextflow.cloud.aws.util

import software.amazon.nio.spi.s3.NextflowS3Path
import software.amazon.nio.spi.s3.S3Path
import software.amazon.nio.spi.s3.S3PathFactory
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class S3PathFactoryTest extends Specification {

    def 'should parse s3 paths' () {

        when:
        def path = S3PathFactory.parse(S3_PATH)
        then:
        S3PathFactory.isS3Path(path)
        with(path as NextflowS3Path) {
            toS3Path().bucketName() == BUCKET
            toS3Path().getKey() == KEY
        }

        when:
        def str = S3PathFactory.getUriString(path)
        then:
        str == S3_PATH

        
        where:
        S3_PATH                                                 | BUCKET        | KEY
        's3://cbcrg-eu/raw/x_r1.fq'                             | 'cbcrg-eu'    | 'raw/x_r1.fq'
        's3://cbcrg-eu/raw/**_R1*{fastq,fq,fastq.gz,fq.gz}'     | 'cbcrg-eu'    | 'raw/**_R1*{fastq,fq,fastq.gz,fq.gz}'

    }

    def 'should ignore double slashes' () {
        when:
        def path = S3PathFactory.parse('s3://cbcrg-eu/raw//x_r1.fq' )
        then:
        S3PathFactory.getUriString(path) == 's3://cbcrg-eu/raw/x_r1.fq'
    }
}
