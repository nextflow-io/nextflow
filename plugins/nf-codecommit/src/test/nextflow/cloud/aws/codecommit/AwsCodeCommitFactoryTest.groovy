/*
 * Copyright 2013-2024, Seqera Labs
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

package nextflow.cloud.aws.codecommit

import nextflow.scm.GitUrl
import nextflow.scm.ProviderConfig
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AwsCodeCommitFactoryTest extends Specification {

    def 'should create a new provider instance' () {
        given:
        def factory = new AwsCodeCommitFactory()
        and:
        def config = Mock(ProviderConfig)

        when:
        def result = factory.createProviderInstance(config, 'some name')
        then:
        result == null
        and:
        config.getPlatform() >> 'any'

    }

    def 'should create config instance' () {
        given:
        def factory = new AwsCodeCommitFactory()

        when:
        def result = factory.createConfigInstance('foo', [:])
        then:
        result == null

        when:
        result = factory.createConfigInstance('codecommit', [:])
        then:
        result instanceof AwsCodeCommitProviderConfig
        result.platform == 'codecommit'
        result.name == 'codecommit'

        when:
        result = factory.createConfigInstance('my-aws-repo', [platform: 'codecommit', user:'foo', password:'xyz'])
        then:
        result instanceof AwsCodeCommitProviderConfig
        result.platform == 'codecommit'
        result.name == 'my-aws-repo'
        result.user == 'foo'
        result.password == 'xyz'

    }

    def 'should get a config' () {
        given:
        def factory = new AwsCodeCommitFactory()

        when:
        def configs = [
                new ProviderConfig('github', [:]),
                new AwsCodeCommitProviderConfig('codecommit', [platform:'codecommit'])
        ]
        and:
        def result = factory.getConfig(configs, new GitUrl('https://github.com/this/that'))
        then:
        result == null

        /*
         * A CodeCommit config is given without any server specification
         * => attributes should be taken and server path updated
         */
        when:
        configs = [
                new ProviderConfig('github', [:]),
                new AwsCodeCommitProviderConfig('codecommit', [platform:'codecommit', 'user':'foo', password: 'xxx'])
        ]
        and:
        result = factory.getConfig(configs, new GitUrl('https://git-codecommit.eu-west-1.amazonaws.com/v1/repos/my-repo'))
        then:
        result instanceof AwsCodeCommitProviderConfig
        and:
        result.name == 'codecommit'
        result.platform == 'codecommit'
        result.server == 'https://git-codecommit.eu-west-1.amazonaws.com'
        result.region == 'eu-west-1'
        result.user == 'foo'
        result.password == 'xxx'
        and:
        // just modifies the instance in the provided list
        result in configs


        /*
         * No config is given => a new one should be created
         */
        when:
        configs = []
        and:
        result = factory.getConfig(configs, new GitUrl('https://git-codecommit.eu-west-1.amazonaws.com/v1/repos/my-repo'))
        then:
        result instanceof AwsCodeCommitProviderConfig
        and:
        result.name == 'codecommit'
        result.platform == 'codecommit'
        result.server == 'https://git-codecommit.eu-west-1.amazonaws.com'
        result.region == 'eu-west-1'
        and:
        // creates a new instance
        result !in configs

        /*
         * A CodeCommit config with a matching server is given => it should be used
         */
        when:
        configs = [
                new ProviderConfig('github', [:]),
                new AwsCodeCommitProviderConfig('my-repo', [platform:'codecommit', server: 'https://git-codecommit.eu-west-1.amazonaws.com'])
        ]
        and:
        result = factory.getConfig(configs, new GitUrl('https://git-codecommit.eu-west-1.amazonaws.com/v1/repos/my-repo'))
        then:
        result instanceof AwsCodeCommitProviderConfig
        and:
        result.name == 'my-repo'
        result.platform == 'codecommit'
        result.server == 'https://git-codecommit.eu-west-1.amazonaws.com'
        result.region == 'eu-west-1'
        and:
        // just modifies the instance in the provided list
        result in configs

        /*
         * A CodeCommit config is given for a different region/server
         * => a new config should be created
         */
        when:
        configs = [
                new ProviderConfig('github', [:]),
                new AwsCodeCommitProviderConfig('my-repo', [platform:'codecommit', server: 'https://git-codecommit.us-east-1.amazonaws.com', user: 'foo'])
        ]
        and:
        result = factory.getConfig(configs, new GitUrl('https://myself@git-codecommit.eu-west-1.amazonaws.com/v1/repos/my-repo'))
        then:
        result instanceof AwsCodeCommitProviderConfig
        and:
        result.name == 'codecommit'
        result.platform == 'codecommit'
        result.server == 'https://git-codecommit.eu-west-1.amazonaws.com'
        result.region == 'eu-west-1'
        result.user == 'myself'
        and:
        result !in configs
    }
}
