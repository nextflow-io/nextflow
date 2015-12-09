package nextflow.scm

import nextflow.exception.AbortOperationException

/**
 * Created by mlogghe on 3/11/2015.
 */
class StashRepositoryProvider extends RepositoryProvider{
    protected String projectKey
    protected String repositorySlug

    StashRepositoryProvider(String project, ProviderConfig config=null) {
        this.project = project
        def items = project.split('/')
        this.projectKey = items[0].toUpperCase()
        this.repositorySlug = items[1]

        this.config = config ?: new ProviderConfig('stash')
    }

    @Override
    String getName() {
        "Stash"
    }

    @Override
    String getEndpointUrl() {
        return "${config.endpoint}/rest/api/1.0/projects/${projectKey}/repos/${repositorySlug}"
    }

    String getDefaultBranch() {
        invokeAndParseResponse("${getEndpointUrl()}/branches/default") ?. displayId
    }

    @Override
    String getContentUrl(String path) {
        "${config.server}/projects/${projectKey}/repos/${repositorySlug}/browse/${path}?at=${getDefaultBranch()}&raw"
    }

    @Override
    String getCloneUrl() {
        Map response = invokeAndParseResponse( getEndpointUrl() )

        if (response?.scmId != "git"){
            throw new AbortOperationException("Stash repository at ${getRepositoryUrl()} is not supporting Git")
        }

        def result = response?.links?.clone?.find{ it.name == "http"} as Map
        if( !result )
            throw new IllegalStateException("Missing clone URL for: $project")

        return result.href
    }

    @Override
    String getRepositoryUrl() {
        "${config.server}/$project"
    }

    @Override
    protected byte[] readBytes(String path) {
        def url = getContentUrl(path)
        def response  = invoke(url)
        response
    }
}
