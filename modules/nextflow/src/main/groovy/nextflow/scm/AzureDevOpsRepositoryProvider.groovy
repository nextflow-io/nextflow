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
 */

package nextflow.scm

import groovy.transform.CompileDynamic
import groovy.transform.CompileStatic
import groovy.transform.Memoized

/**
 * Implements a repository provider for Azure DevOps service
 *
 * @author Tobias Neumann <tobias.neumann.at@gmail.com>
 */
@CompileStatic
final class AzureDevOpsRepositoryProvider extends RepositoryProvider {

    AzureDevOpsRepositoryProvider(String project, ProviderConfig config=null) {
        this.project = project
        this.config = config ?: new ProviderConfig('azuredevops')
    }

    /** {@inheritDoc} */
    @Override
    String getName() { "Azure DevOps" }

    /** {@inheritDoc} */
    @Override
    String getEndpointUrl() {
        //https://dev.azure.com/quantrocode/slamseq/_apis/git/repositories/slamseq
        "https://dev.azure.com/quantrocode/slamseq/_apis/git/repositories/slamseq"
        //"${config.endpoint}/${project}/_apis/projects?api-version=2.0"
    }

    /** {@inheritDoc} */
    @Override
    String getContentUrl( String path ) {
        "https://dev.azure.com/quantrocode/slamseq/_apis/git/repositories/slamseq/items?download=false&includeContent=true&includeContentMetadata=false&api-version=6.0&$format=json&path=$path"
        //"${config.endpoint}/repos/$project/contents/$path"
    }

    /** {@inheritDoc} */
    @Override
    String getCloneUrl() {

        //return "https://${organization}.visualstudio.com/${repo}/_git/${repo}"
        return "https://dev.azure.com/quantrocode/slamseq/_git/slamseq"
        /*Map response = invokeAndParseResponse( getEndpointUrl() )

        def result = response.get('clone_url')
        if( !result )
            throw new IllegalStateException("Missing clone URL for: $project")

        return result/*
    }

    @Override
    @CompileDynamic
    @Memoized
    List<BranchInfo> getBranches() {
        // https://dev.azure.com/quantrocode/slamseq/_apis/git/repositories/slamseq/refs?api-version=6.0
        final url = "https://dev.azure.com/quantrocode/slamseq/_apis/git/repositories/slamseq/refs?api-version=6.0"
        //final url = "${config.endpoint}/repos/$project/branches"
        this.<BranchInfo>invokeAndResponseWithPaging(url, { Map branch -> new BranchInfo(branch.name as String, branch.commit?.sha as String) })
    }

    @Override
    @CompileDynamic
    @Memoized
    List<TagInfo> getTags() {
        // https://dev.azure.com/quantrocode/slamseq/_apis/git/repositories/slamseq/refs?filter=tags/&api-version=6.0
        final url = "https://dev.azure.com/quantrocode/slamseq/_apis/git/repositories/slamseq/refs?filter=tags/&api-version=6.0"
        // final url = "${config.endpoint}/repos/$project/tags"
        this.<TagInfo>invokeAndResponseWithPaging(url, { Map tag -> new TagInfo(tag.name as String, tag.commit?.sha as String)})
    }

    /** {@inheritDoc} */
    @Override
    String getRepositoryUrl() {
        //"${config.server}/$project"
        "https://dev.azure.com/quantrocode/_git/slamseq"
    }

    /** {@inheritDoc} */
    @Override
    byte[] readBytes(String path) {

        def url = getContentUrl(path)
        Map response  = invokeAndParseResponse(url)
        response.get('content')?.toString()?.decodeBase64()

    }

}
