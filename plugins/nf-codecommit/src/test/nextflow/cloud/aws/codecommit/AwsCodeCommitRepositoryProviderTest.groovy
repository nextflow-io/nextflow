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

import nextflow.scm.RepositoryProvider
import spock.lang.IgnoreIf
import spock.lang.Requires
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@IgnoreIf({System.getenv('NXF_SMOKE')})
@Requires({System.getenv('AWS_ACCESS_KEY_ID') && System.getenv('AWS_SECRET_ACCESS_KEY')})
class AwsCodeCommitRepositoryProviderTest extends Specification {

    def 'should get repo url' () {
        given:
        def config = new AwsCodeCommitProviderConfig('git-codecommit.eu-west-1.amazonaws.com')
        and:
        def provider = new AwsCodeCommitRepositoryProvider('codecommit-eu-west-1/my-repo', config)

        expect:
        provider.getCloneUrl() == 'https://git-codecommit.eu-west-1.amazonaws.com/v1/repos/my-repo'
        and:
        provider.getRepositoryUrl() == "https://git-codecommit.eu-west-1.amazonaws.com/v1/repos/my-repo"
    }

    def 'should read content' () {
        given:
        def config = new AwsCodeCommitProviderConfig('git-codecommit.eu-west-1.amazonaws.com')
        and:
        def provider = new AwsCodeCommitRepositoryProvider('codecommit-eu-west-1/my-repo', config)

        expect:
        provider.readText('main.nf') == '''\
                nextflow.enable.dsl=2
                
                workflow {
                  sayHello()
                }
                
                process sayHello {
                  /echo Hello world/
                }
                '''.stripIndent().rightTrim()
    }

    def 'should read content with revision' () {
        given:
        def config = new AwsCodeCommitProviderConfig('git-codecommit.eu-west-1.amazonaws.com')
        and:
        def provider = new AwsCodeCommitRepositoryProvider('codecommit-eu-west-1/my-repo', config)
        and:
        provider.revision = 'dev1'
        expect:
        provider.readText('main.nf') == 'println "Hello world in dev branch!"\n'
    }

    def 'should fetch repo tags'() {
        given:
        def config = new AwsCodeCommitProviderConfig('git-codecommit.eu-west-1.amazonaws.com')
        and:
        def provider = new AwsCodeCommitRepositoryProvider('codecommit-eu-west-1/my-repo', config)

        when:
        // uses repo at
        // https://git-codecommit.eu-west-1.amazonaws.com/v1/repos/my-repo
        def result = provider.getTags() as Set
        then:
        result == [
                new RepositoryProvider.TagInfo('v0.1', '1a88516b5e382d0d68bfa01c18eab6c2067c0595'),
                new RepositoryProvider.TagInfo('v0.2', 'c673d3d55be190c54db2056690b71e285fe5b3d8')] as Set

    }

    def 'should fetch repo branches'() {
        given:
        def config = new AwsCodeCommitProviderConfig('git-codecommit.eu-west-1.amazonaws.com')
        and:
        def provider = new AwsCodeCommitRepositoryProvider('codecommit-eu-west-1/my-repo', config)

        when:
        // uses repo at
        // https://git-codecommit.eu-west-1.amazonaws.com/v1/repos/my-repo
        def result = provider.getBranches() as Set
        then:
        result == [
                new RepositoryProvider.BranchInfo('master', 'c820e0904d9ce4404e005e3cc910502300b36ba3'),
                new RepositoryProvider.BranchInfo('dev1', 'c90422a1b4823f1c0980bbf8cab261e45a351622')] as Set

    }

}
