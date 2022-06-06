/*
 * Copyright 2020, Seqera Labs
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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


import spock.lang.IgnoreIf
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@IgnoreIf({System.getenv('NXF_SMOKE')})
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
        provider.readText('main.nf') == 'println "Hello world!"\n'
    }

}
