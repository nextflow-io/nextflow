package nextflow.util


import test.BaseSpec
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class SecretHelperTest extends BaseSpec {

    def 'should remove secret key' () {
        expect:
        SecretHelper.secureEnvString('a=b') == 'a=b'
        and:
        SecretHelper.secureEnvString('aws_key=12345') == 'aws_key=[secure]'
        SecretHelper.secureEnvString('AWS_KEY=12345') == 'AWS_KEY=[secure]'

        SecretHelper.secureEnvString('''\
                foo=hello    
                aws_key=d7sds89
                git_token=909s-ds-'''
                .stripIndent() ) ==
                '''\
                foo=hello    
                aws_key=[secure]
                git_token=[secure]'''.stripIndent()

    }

    def 'should remove secrets' () {
        given:
        def obj = [foo: 'hello',
                   awsKey: '1234',
                   githubToken: 'abc' ]

        expect:
        SecretHelper.hideSecrets(obj) == [
                foo: 'hello',
                awsKey: '[secret]',
                githubToken: '[secret]' ]
    }

    def 'should remove nested secrets' () {
        given:
        def obj = [
                aws: [secretKey: 'abc', accessKey: 'zzz'],
                github: [ [token: 'xxx'], [endpoint: 'this is good'] ]
                ]

        expect:
        SecretHelper.hideSecrets(obj) == [
                aws: [accessKey: '[secret]', secretKey: '[secret]'],
                github: [ [token: '[secret]'], [endpoint: 'this is good']
                ]
        ]
    }

}
