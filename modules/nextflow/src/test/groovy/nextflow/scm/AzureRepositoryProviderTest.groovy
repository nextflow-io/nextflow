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

import nextflow.exception.AbortOperationException
import spock.lang.IgnoreIf
import spock.lang.Requires
import spock.lang.Specification
/**
 *
 * @author Tobias Neumann <tobias.neumann.at@gmail.com>
 */
class AzureRepositoryProviderTest extends Specification {

    static final String CONFIG = '''
        providers {

            azurerepos {
                server = 'https://dev.azure.com'
                endpoint = 'https://dev.azure.com'
                platform = 'azurerepos'
                user = 'myname'
                password = 'mypassword'
                token = '65eaa9c8ef52460d22a93307fe0aee76289dc675'
                revision = 'master'
            }

        }
        '''

    def 'should parse repo fields from path' () {
        expect:
        AzureRepositoryProvider.getUniformPath(PATH) == EXPECTED

        where:
        PATH                                | EXPECTED
        't-neumann/hello'                   | ['t-neumann', 'hello', 'hello']
        'ORGANIZATION/PROJECT/hello'        | ['ORGANIZATION','PROJECT','hello']
        'ORGANIZATION/PROJECT/_git/hello'   | ['ORGANIZATION','PROJECT','hello']

    }

    def 'should throw exception if wrong path' () {
        when:
        def path = AzureRepositoryProvider.getUniformPath(PATH)

        then :
        def exception = thrown(IllegalArgumentException)
        exception?.message == EXCEPTION

        where:
        PATH                                | EXCEPTION
        'incorrect_path_1'                  | "Unexpected Azure repository path format - offending value: 'incorrect_path_1'"
        'ORG/PROJ/hello/incorrect'          | "Unexpected Azure repository path format - offending value: 'ORG/PROJ/hello/incorrect'"
    }

    def 'should return repo url' () {

        given:
        def config = new ConfigSlurper().parse(CONFIG)
        def obj = new ProviderConfig('azurerepos', config.providers.azurerepos as ConfigObject)

        expect:
        new AzureRepositoryProvider(PATH, obj).getEndpointUrl() == EXPECTED

        where:
        PATH                                | EXPECTED
        't-neumann/hello'                   | 'https://dev.azure.com/t-neumann/hello/_apis/git/repositories/hello'
        'ORGANIZATION/PROJECT/hello'        | 'https://dev.azure.com/ORGANIZATION/PROJECT/_apis/git/repositories/hello'
        'ORGANIZATION/PROJECT/_git/hello'   | 'https://dev.azure.com/ORGANIZATION/PROJECT/_apis/git/repositories/hello'
    }

    def 'should return project URL' () {

        given:
        def config = new ConfigSlurper().parse(CONFIG)
        def obj = new ProviderConfig('azurerepos', config.providers.azurerepos as ConfigObject)

        expect:
        new AzureRepositoryProvider(PATH, obj).getRepositoryUrl() == EXPECTED
        where:
        PATH                                | EXPECTED
        't-neumann/hello'                   | 'https://dev.azure.com/t-neumann/hello'
        'ORGANIZATION/PROJECT/hello'        | 'https://dev.azure.com/ORGANIZATION/PROJECT/hello'
        'ORGANIZATION/PROJECT/_git/hello'   | 'https://dev.azure.com/ORGANIZATION/PROJECT/_git/hello'

    }

    def 'should return content URL' () {

        given:
        def config = new ConfigSlurper().parse(CONFIG)
        def obj = new ProviderConfig('azurerepos', config.providers.azurerepos as ConfigObject)

        expect:
        new AzureRepositoryProvider('t-neumann/hello', obj).getContentUrl('main.nf') == 'https://dev.azure.com/t-neumann/hello/_apis/git/repositories/hello/items?download=false&includeContent=true&includeContentMetadata=false&api-version=6.0&\$format=json&path=main.nf'

    }

    def 'should return content URL for revision' () {

        given:
        def config = new ConfigSlurper().parse(CONFIG)
        def obj = new ProviderConfig('azurerepos', config.providers.azurerepos as ConfigObject)

        expect:
        new AzureRepositoryProvider('t-neumann/hello', obj).setRevision("a-branch").getContentUrl('main.nf') == 'https://dev.azure.com/t-neumann/hello/_apis/git/repositories/hello/items?download=false&includeContent=true&includeContentMetadata=false&api-version=6.0&\$format=json&path=main.nf&versionDescriptor.version=a-branch'

    }

    /*
     * The NXF_AZURE_REPOS_TOKEN must contain `user_name : password` provided by the Azure repos service.
     *
     * Check for *Generate Git credentials* in the *Clone Repository* box
     */
    @IgnoreIf({System.getenv('NXF_SMOKE')})
    @Requires({System.getenv('NXF_AZURE_REPOS_TOKEN')})
    def 'should read file content'() {
        given:
        def token = System.getenv('NXF_AZURE_REPOS_TOKEN')
        def config = new ProviderConfig('azurerepos').setAuth(token)

        when:
        // uses repo at
        //  https://pditommaso.visualstudio.com/nf-azure-repo/_git/nf-azure-repo
        def repo = new AzureRepositoryProvider('pditommaso/nf-azure-repo', config)
        def result = repo.readText('main.nf')
        then:
        result == 'println "Hello from Azure repos!"'

    }

    @IgnoreIf({System.getenv('NXF_SMOKE')})
    @Requires({System.getenv('NXF_AZURE_REPOS_TOKEN')})
    def 'should fetch repo tags'() {
        given:
        def token = System.getenv('NXF_AZURE_REPOS_TOKEN')
        def config = new ProviderConfig('azurerepos').setAuth(token)

        when:
        // uses repo at
        //  https://pditommaso.visualstudio.com/nf-azure-repo/_git/nf-azure-repo
        def repo = new AzureRepositoryProvider('pditommaso/nf-azure-repo', config)
        def result = repo.getTags()
        then:
        result == [
                new RepositoryProvider.TagInfo('v1.0', '769afa02858ef904756eb5b9a9a04ce8ef869f42'),
                new RepositoryProvider.TagInfo('v2.0', '1f73719675872801cc32fa7fc4ac7f2c04320714')]

    }

    @IgnoreIf({System.getenv('NXF_SMOKE')})
    @Requires({System.getenv('NXF_AZURE_REPOS_TOKEN')})
    def 'should fetch repo branches'() {
        given:
        def token = System.getenv('NXF_AZURE_REPOS_TOKEN')
        def config = new ProviderConfig('azurerepos').setAuth(token)

        when:
        // uses repo at
        //  https://pditommaso.visualstudio.com/nf-azure-repo/_git/nf-azure-repo
        def repo = new AzureRepositoryProvider('pditommaso/nf-azure-repo', config)
        def result = repo.getBranches()
        then:
        result == [
                new RepositoryProvider.BranchInfo('dev', 'cc0ca18640a5c995231e22d91f1527d5155d024b'),
                new RepositoryProvider.BranchInfo('feature-x', '13456a001ba5a27d643755614ab8e814d94ef888'),
                new RepositoryProvider.BranchInfo('master', 'f84130388714582e20f0e2ff9a44b41978ec8929'),
        ]
    }

    @IgnoreIf({System.getenv('NXF_SMOKE')})
    @Requires({System.getenv('NXF_AZURE_REPOS_TOKEN')})
    def 'should not fetch repo file from non existing revision'() {
        given:
        def token = System.getenv('NXF_AZURE_REPOS_TOKEN')
        def config = new ProviderConfig('azurerepos').setAuth(token)

        when:
        // uses repo at
        //  https://pditommaso.visualstudio.com/nf-azure-repo/_git/nf-azure-repo
        def repo = new AzureRepositoryProvider('pditommaso/nf-azure-repo', config)
        repo.revision = 'a_fake_branch'
        def result = repo.readText('file1.txt')

        then:
        thrown(AbortOperationException)
    }

    @IgnoreIf({System.getenv('NXF_SMOKE')})
    @Requires({System.getenv('NXF_AZURE_REPOS_TOKEN')})
    def 'should fetch repo file from revision'() {
        given:
        def token = System.getenv('NXF_AZURE_REPOS_TOKEN')
        def config = new ProviderConfig('azurerepos').setAuth(token)

        when:
        // uses repo at
        //  https://pditommaso.visualstudio.com/nf-azure-repo/_git/nf-azure-repo
        def repo = new AzureRepositoryProvider('pditommaso/nf-azure-repo', config)
        repo.revision = 'dev'
        def result = repo.readText('file-on-dev.txt')
        then:
        result=='hello\n'
    }

    @IgnoreIf({System.getenv('NXF_SMOKE')})
    @Requires({System.getenv('NXF_AZURE_REPOS_TOKEN')})
    def 'should infer revision type for commit revision'() {
        given:
        def token = System.getenv('NXF_AZURE_REPOS_TOKEN')
        def config = new ProviderConfig('azurerepos').setAuth(token)

        when:
        // uses repo at
        //  https://pditommaso.visualstudio.com/nf-azure-repo/_git/nf-azure-repo
        def repo = new AzureRepositoryProvider('pditommaso/nf-azure-repo', config)
        repo.revision = 'cc0ca18640a5c995231e22d91f1527d5155d024b'
        def result = repo.readText('file-on-dev.txt')
        then:
        result=='hello\n'
    }
}
