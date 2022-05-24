/*
 * Copyright 2020-2022, Seqera Labs
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
 * Implements a repository provider for Azure Repos service
 *
 * See https://azure.microsoft.com/en-us/services/devops/repos/
 *
 * @author Tobias Neumann <tobias.neumann.at@gmail.com>
 */
@CompileStatic
final class AzureRepositoryProvider extends RepositoryProvider {

    private String user
    private String repo
    private String continuationToken

    AzureRepositoryProvider(String project, ProviderConfig config=null) {
        /*
        Azure repo format follows Organization/Project/Repository where Project can be optional
        If Project is not present then Repository is used as Project (and also as Repository)
         */
        def tokens = project.tokenize('/')
        this.repo = tokens.removeLast()
        if( tokens.size() == 1){
            this.project = [tokens.first(), this.repo].join('/')
        }else{
            this.project = tokens.join('/')
        }
        this.config = config ?: new ProviderConfig('azurerepos')
        this.continuationToken = null
    }

    /** {@inheritDoc} */
    @Override
    String getName() { "Azure Repos" }

    /** {@inheritDoc} */
    @Override
    String getEndpointUrl() {
        "${config.endpoint}/${project}/_apis/git/repositories/${repo}"
    }

    /** {@inheritDoc} */
    @Override
    String getContentUrl( String path ) {
        // see
        // https://docs.microsoft.com/en-us/rest/api/azure/devops/git/items/get?view=azure-devops-rest-6.0
        //
        def queryParams =[
                'download':false,
                'includeContent':true,
                'includeContentMetadata':false,
                "api-version":6.0,
                '$format':'json',
                'path':path
        ] as Map<String,Object>
        if( revision )
            queryParams['versionDescriptor.version']=revision
        def queryString = queryParams.collect({ "$it.key=$it.value"}).join('&')
        def result = "$endpointUrl/items?$queryString"
        result
    }

    @Override
    @CompileDynamic
    List<BranchInfo> getBranches() {
        // https://docs.microsoft.com/en-us/rest/api/azure/devops/git/refs/list?view=azure-devops-rest-6.0#refs-heads
        final url = "$endpointUrl/refs?filter=heads&api-version=6.0"
        this.<BranchInfo>invokeAndResponseWithPaging(url, { Map branch ->
            new BranchInfo(strip(branch.name), branch.objectId as String)
        })
    }

    @Override
    @CompileDynamic
    @Memoized
    List<TagInfo> getTags() {
        // https://docs.microsoft.com/en-us/rest/api/azure/devops/git/refs/list?view=azure-devops-rest-6.0#refs-tags
        final url = "$endpointUrl/refs?filter=tags&api-version=6.0"
        this.<TagInfo>invokeAndResponseWithPaging(url, { Map tag ->
            new TagInfo(strip(tag.name), tag.objectId as String)
        })
    }

    private String strip(value) {
        if( !value )
            return value
        return value
                .toString()
                .replace('refs/tags/','')
                .replace('refs/heads/','')
    }

    @Memoized
    protected <T> List<T> invokeAndResponseWithPaging(String url, Closure<T> parse) {
        final result = new ArrayList<T>(50)
        do {
            final resp = invokeAndParseResponse(url)
            final tags = resp.value as List<Map>
            for( Map entry : tags ) {
                final item = parse.call(entry)
                if( item )
                    result.add(item)
            }

            // the continuation token is extract when parsing the response
            // see the method `checkResponse`
            url = continuationToken  ? "${url}&continuationToken=${continuationToken}" : null
        }
        while( url != null )
        return result
    }

    /**
     * Check for response error status. Throws a {@link nextflow.exception.AbortOperationException} exception
     * when a 401 or 403 error status is returned.
     *
     * @param connection A {@link HttpURLConnection} connection instance
     */
    protected checkResponse( HttpURLConnection connection ) {

        if (connection.getHeaderFields().containsKey("x-ms-continuationtoken")) {
            this.continuationToken = connection.getHeaderField("x-ms-continuationtoken");
        } else {
            this.continuationToken = null
        }

        super.checkResponse(connection)
    }

    /** {@inheritDoc} */
    @Override
    String getCloneUrl() {
        return "https://dev.azure.com/${project}/_git/${repo}"
    }

    /** {@inheritDoc} */
    @Override
    String getRepositoryUrl() {
        "${config.server}/$project"
    }

    /** {@inheritDoc} */
    @Override
    byte[] readBytes(String path) {
        final url = getContentUrl(path)
        final response = invokeAndParseResponse(url)
        return response.get('content')?.toString()?.getBytes()
    }

}
