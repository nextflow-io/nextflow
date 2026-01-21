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
 */

package nextflow.scm

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DefaultRepositoryFactoryLoaderTest extends Specification {

    def 'should create repository provider object' () {

        def provider

        when:
        provider = DefaultRepositoryFactoryLoader.instance.createRepositoryProvider(new ProviderConfig('github'),'project/x')
        then:
        provider instanceof GithubRepositoryProvider
        provider.endpointUrl == 'https://api.github.com/repos/project/x'

        when:
        provider = DefaultRepositoryFactoryLoader.instance.createRepositoryProvider(new ProviderConfig('gitlab'),'project/y')
        then:
        provider instanceof GitlabRepositoryProvider
        provider.endpointUrl == 'https://gitlab.com/api/v4/projects/project%2Fy'

        when:
        provider = DefaultRepositoryFactoryLoader.instance.createRepositoryProvider(new ProviderConfig('bitbucket'),'project/z')
        then:
        provider instanceof BitbucketRepositoryProvider
        provider.endpointUrl == 'https://api.bitbucket.org/2.0/repositories/project/z'

        when:
        provider = DefaultRepositoryFactoryLoader.instance.createRepositoryProvider(new ProviderConfig('local', [path:'/user/data']),'local/w')
        then:
        provider.endpointUrl == 'file:/user/data/w'
    }

}
