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

package io.nextflow.gradle.tasks

import com.google.gson.Gson
import groovy.transform.CompileStatic
import groovy.transform.Memoized
import io.nextflow.gradle.util.GithubClient
import org.apache.commons.codec.digest.DigestUtils
import org.gradle.api.DefaultTask
import org.gradle.api.GradleException
import org.gradle.api.provider.ListProperty
import org.gradle.api.provider.Property
import org.gradle.api.tasks.Input
import org.gradle.api.tasks.TaskAction
/**
 * Upload project artifacts to the corresponding Github repository releases page
 *
 * Based on https://github.com/mgk/s3-plugin/blob/master/src/main/groovy/com/github/mgk/gradle/S3Plugin.groovy
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class GithubUploader extends DefaultTask {

    /**
     * The source file to upload
     */
    @Input final ListProperty<String> assets = project.objects.listProperty(String)

    @Input final Property<String> repo = project.objects.property(String)

    @Input final Property<String> release = project.objects.property(String)

    @Input final Property<Boolean> unstable = project.objects.property(Boolean)

    @Input String owner

    @Input boolean overwrite = false

    @Input boolean dryRun = false

    @Input boolean skipExisting

    @Input boolean draft

    @Input String authToken

    @Input String userName

    @Input boolean ignore

    @Memoized
    private GithubClient getClient() {
        new GithubClient(authToken: authToken, owner: owner, repo: repo.get(), userName: userName)
    }

    @TaskAction
    def task() {
        final files = (List) assets.get()
        for( String it : files ) {
            upload( new File(it) )
        }
    }

    private void upload(File sourceFile) {
        if( ignore )
            return

        if( !sourceFile.exists() )
            throw new GradleException("Github upload failed -- source file does not exist: $sourceFile")

        final fileName = sourceFile.name
        final asset = client.getReleaseAsset(release.get(), fileName)
        if ( asset ) {
            if( skipExisting && isSame(sourceFile, asset) ) {
                logger.quiet("${owner}/${repo.get()}/${fileName} already exists -- skipping")
            }
            else if (overwrite) {
                updateRelease(sourceFile)
            }
            else {
                throw new GradleException("${owner}/${repo.get()}/${fileName} already exists -- overwrite refused")
            }
        }
        else {
            uploadRelease(sourceFile)
        }
    }

    private void updateRelease(File sourceFile) {
        if( dryRun ) {
            logger.quiet("Will update ${sourceFile} → github.com://$owner/${repo.get()}")
        }
        else {
            logger.quiet("Updating ${sourceFile} → github.com://$owner/${repo.get()}")
            client.deleteReleaseAsset(release.get(), sourceFile.name)
            client.uploadReleaseAsset(release.get(), sourceFile, mime(sourceFile.name))
        }
    }

    private void uploadRelease(File sourceFile) {
        if( dryRun ) {
            logger.quiet("Will upload ${sourceFile} → github.com://$owner/${repo.get()}")
        }
        else {
            logger.quiet("Uploading ${sourceFile} → github.com://$owner/${repo.get()}")
            final rel = client.getRelease(release.get()) ?: client.createRelease(release.get(), unstable.get())
            final releaseId = (rel.id as Long).toString()
            client.uploadReleaseAsset(releaseId, sourceFile, mime(sourceFile.name))
        }
    }

    private String mime(String fileName) {
        if( fileName.endsWith('.zip') ) {
            return "application/zip"
        }
        if( fileName.endsWith('.json') ) {
            return "application/json"
        }
        throw new IllegalArgumentException("Unknown file type: $fileName")
    }

    private boolean isSame(File sourceFile, InputStream asset ) {
        sourceFile.name.endsWith('.json')
            ? isSameJson0(sourceFile, asset)
            : isSameBin0(sourceFile, asset)
    }

    private boolean isSameBin0(File sourceFile, InputStream asset ) {
        final d1 = sourceFile
                .withInputStream { InputStream it -> DigestUtils.sha512Hex(it) }
        final d2 = DigestUtils.sha512Hex(asset)
        return d1 == d2
    }

    private boolean isSameJson0(File sourceFile, InputStream asset) {
        def gson = new Gson()
        def j1 = gson.fromJson(sourceFile.text, Map)
        def j2 = gson.fromJson( new InputStreamReader(asset), Map)
        if( j1.version != j2.version ) {
            logger.quiet("Plugin metafile $sourceFile does not match versions: local=$j1.version; remote: $j2.version")
            return false
        }
        if( j1.url != j2.url ) {
            logger.quiet("Plugin metafile $sourceFile does not match urls: local=$j1.url; remote: $j2.url")
            return false
        }
        if( j1.sha512sum != j2.sha512sum ) {
            logger.quiet("Plugin metafile $sourceFile does not match sha512sum: local=$j1.sha512sum; remote: $j2.sha512sum")
            return false
        }
        return true
    }

}
