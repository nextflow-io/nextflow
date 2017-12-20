package nextflow.executor

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AwsOptionsTest extends Specification {

    def 'should return aws cli' () {

        given:
        AwsOptions opts

        when:
        opts = new AwsOptions()
        then:
        opts.awsCli == 'aws'

        when:
        opts = new AwsOptions(cliPath: '/foo/bin/aws')
        then:
        opts.awsCli == '/foo/bin/aws'

        when:
        opts = new AwsOptions(cliPath: '/foo/bin/aws', region: 'eu-west-1')
        then:
        opts.awsCli == '/foo/bin/aws --region eu-west-1'
    }
}
