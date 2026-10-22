package nextflow.util


import spock.lang.Timeout
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

    def 'should hide secrets with the same name rule as stripSecrets' () {
        given:
        def obj = [
                azure: [storage: [accountName: 'store', accountKey: 'k1'], batch: [accountKey: 'k2']],
                tower: [license: 'l1', basicAuth: 'a1'],
                db: [passwd: 'p1', credentials: 'c1', privateKey: 'pk'],
                aws: [region: 'eu-west-1'] ]

        expect:
        SecretHelper.hideSecrets(obj) == [
                azure: [storage: [accountName: 'store', accountKey: '[secret]'], batch: [accountKey: '[secret]']],
                tower: [license: '[secret]', basicAuth: '[secret]'],
                db: [passwd: '[secret]', credentials: '[secret]', privateKey: '[secret]'],
                aws: [region: 'eu-west-1'] ]
    }

    def 'should secure env lines with the same name rule as stripSecrets' () {
        expect:
        SecretHelper.secureEnvString('DB_PASSWORD=p1') == 'DB_PASSWORD=[secure]'
        SecretHelper.secureEnvString('AZURE_STORAGE_ACCOUNT_KEY=k1==') == 'AZURE_STORAGE_ACCOUNT_KEY=[secure]'
        SecretHelper.secureEnvString('TOWER_LICENSE=l1') == 'TOWER_LICENSE=[secure]'
        SecretHelper.secureEnvString('GITHUB_CREDENTIALS=c1') == 'GITHUB_CREDENTIALS=[secure]'
        and: 'the whole value is masked, also when it contains `=`'
        SecretHelper.secureEnvString('AZURE_STORAGE_SAS_TOKEN=sv=2020&sig=abc') == 'AZURE_STORAGE_SAS_TOKEN=[secure]'
        and: 'a secret after another assignment on the same line'
        SecretHelper.secureEnvString('FOO=bar MY_TOKEN=xyz') == 'FOO=bar MY_TOKEN=[secure]'
        and: 'not secrets'
        SecretHelper.secureEnvString('FOO=bar\nKEY_ID=x\nREGION=eu') == 'FOO=bar\nKEY_ID=x\nREGION=eu'
    }

    def 'should not hide manifest values' () {
        given:
        def obj = [manifest: [author: 'Paolo', license: 'MIT', name: 'foo'], tower: [license: 'l1']]

        expect:
        SecretHelper.hideSecrets(obj) == [manifest: [author: 'Paolo', license: 'MIT', name: 'foo'], tower: [license: '[secret]']]
    }

    def 'should secure env lines with any name format' () {
        expect:
        SecretHelper.secureEnvString('export GH_TOKEN=x') == 'export GH_TOKEN=[secure]'
        SecretHelper.secureEnvString('  NXF_TOKEN = x') == '  NXF_TOKEN =[secure]'
        SecretHelper.secureEnvString('aws.region=eu-west-1') == 'aws.region=[secure]'
        SecretHelper.secureEnvString('AWS=x') == 'AWS=[secure]'
        SecretHelper.secureEnvString('MY_TOKEN=a FOO=b\nBAR=c') == 'MY_TOKEN=[secure]\nBAR=c'
        SecretHelper.secureEnvString('A=1\nB_SECRET=2\nC=3\nD_PASSWORD=4') == 'A=1\nB_SECRET=[secure]\nC=3\nD_PASSWORD=[secure]'
        SecretHelper.secureEnvString('') == ''
    }

    @Timeout(5)
    def 'should secure env strings in linear time' () {
        expect:
        [ 'FOO=' + 'a' * 100_000, 'FOO=' + 'token' * 20_000, 'a-' * 50_000, '=' * 100_000 ].every { SecretHelper.secureEnvString(it) == it }
    }

}
