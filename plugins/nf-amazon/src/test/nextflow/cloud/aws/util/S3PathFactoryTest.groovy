package nextflow.cloud.aws.util

import com.upplication.s3fs.S3Path
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
        path instanceof S3Path
        with(path as S3Path) {
            getBucket() == BUCKET
            getKey() == KEY
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
