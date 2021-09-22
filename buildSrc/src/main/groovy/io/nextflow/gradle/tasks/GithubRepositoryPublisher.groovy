package io.nextflow.gradle.tasks


import com.google.gson.Gson
import com.google.gson.GsonBuilder
import com.google.gson.reflect.TypeToken
import groovy.transform.CompileStatic
import io.nextflow.gradle.model.PluginMeta
import io.nextflow.gradle.model.PluginRelease
import io.nextflow.gradle.util.GithubClient
import org.gradle.api.DefaultTask
import org.gradle.api.GradleException
import org.gradle.api.tasks.Input
import org.gradle.api.tasks.Optional
import org.gradle.api.tasks.TaskAction
/**
 * This task traverse a S3 plugins repo and creates
 * and updates plugins repository index. Finally push
 * the updated index to the Github repository.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class GithubRepositoryPublisher extends DefaultTask {

    @Input List<String> repos

    /**
     * The target plugins repository index HTTP URL
     */
    @Input String indexUrl

    /**
     * The auth access token to post to access Github plugins repo
     */
    @Input String githubToken

    /**
     * The Github user name posting the changes
     */
    @Input String githubUser

    /**
     * The Github email signing the index update
     */
    @Input String githubEmail

    /**
     * The AWS creds profile
     */
    @Input @Optional String profile

    /**
     * The AWS plugins repo region
     */
    @Input @Optional String region

    @Input @Optional Boolean overwrite

    @Input String owner


    String mergeIndex(List<PluginMeta> mainIndex, Map<String,List<PluginRelease>> pluginsToPublish) {

        for( Map.Entry<String,List<PluginRelease>> item : pluginsToPublish ) {
            final pluginId = item.key
            final pluginReleases = item.value
            final PluginMeta indexEntry = mainIndex.find { PluginMeta meta -> meta.id == pluginId }

            if (!indexEntry) {
                mainIndex.add(new PluginMeta(id: pluginId, releases: pluginReleases))
            }
            else {
                for (PluginRelease rel : pluginReleases) {
                    // check if this version already exist in the index
                    final index = indexEntry.releases.findIndexOf { PluginRelease it -> it.version == rel.version }
                    final indexRel = index!=-1 ? indexEntry.releases[index] : null as PluginRelease

                    // if not exists, add to the index
                    if( !indexRel ) {
                        indexEntry.releases << rel
                    }
                    // otherwise verify the checksum matches
                    else if( indexRel.sha512sum != rel.sha512sum ) {
                        if( overwrite ) {
                            indexEntry.releases[index] = rel
                        }
                        else {
                            def msg = "Plugin $pluginId@${rel.version} invalid checksum:\n"
                            msg += "- index sha512sum: $indexRel.sha512sum\n"
                            msg += "- repo sha512sum : $rel.sha512sum\n"
                            msg += "- repo url       : $rel.url"
                            throw new GradleException(msg)
                        }
                    }
                }
            }
        }

        new GsonBuilder()
                .setPrettyPrinting()
                .disableHtmlEscaping()
                .create()
                .toJson(mainIndex)
    }

    List<PluginMeta> parseMainIndex(GithubClient github, String path) {
        // get main repo index
        final indexJson = github.getContent(path)
        final type = new TypeToken<ArrayList<PluginMeta>>(){}.getType()
        return new Gson().fromJson(indexJson, type)
    }


    /*
     * Traverse a S3 bucket and return a map given all releases for each
     * plugin id
     *
     * @return The map holding the plugin releases for each plugin id
     */
    Map<String,List<PluginRelease>> listPlugins() {
        Map<String,List<PluginRelease>> result = [:]
        for( String it : repos ) {
            final rel = getReleases(it)
            if( rel ) result[it] = [rel]
        }
        return result
    }

    PluginRelease getReleases(String repo) {

        final client = new GithubClient(
                authToken: githubToken,
                userName: githubUser,
                owner: owner,
                repo: repo)

        // https://docs.github.com/en/free-pro-team@latest/rest/reference/repos#get-the-latest-release
        final resp = client.latestRelease()
        if( !resp ) {
            logger.quiet("WARN: No release found for repo $owner/$repo")
            return null
        }

        final version = resp.tag_name as String
        if( !version ) {
            logger.quiet("WARN: No version found for repo $owner/$repo")
            return null
        }

        logger.quiet("Merging $owner/$repo@$version")
        final metaFile = "${repo}-${version}-meta.json"
        final json = client.getReleaseAsset(version, metaFile)?.text
        if( !json )
            throw new GradleException("Can't load plugin release metafile $metaFile")
        return new Gson().fromJson(json, PluginRelease)
    }

    @TaskAction
    def apply() {

        // parse indexUrl 
        final gitUrl = new URL(indexUrl)
        final tokns = gitUrl.path.tokenize('/')
        final githubOrg = tokns[0]
        final githubRepo = tokns[1]
        final githubBranch = tokns[2]
        final targetFileName = tokns[3]

        // init github client
        final github = new GithubClient()
        github.userName = githubUser
        github.authToken = githubToken
        github.branch = githubBranch ?: 'master'
        github.repo = githubRepo
        github.owner = githubOrg
        github.email = githubEmail

        // list plugins in the nextflow s3 releases
        logger.quiet("Fetching plugins $repos")
        final pluginsToPublish = listPlugins()

        // fetch the plugins public index
        logger.quiet("Parsing current index $indexUrl")
        def mainIndex = parseMainIndex(github, targetFileName)

        // merge indexes
        logger.quiet("Merging index")
        final result = mergeIndex(mainIndex, pluginsToPublish)

        // push to github
        logger.quiet("Publish merged index to $indexUrl")

        github.pushChange(targetFileName, result.toString() + '\n', "Nextflow plugins update")
    }

}
