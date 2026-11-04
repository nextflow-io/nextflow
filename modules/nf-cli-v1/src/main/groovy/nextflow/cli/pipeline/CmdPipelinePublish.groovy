/*
 * Copyright 2013-2026, Seqera Labs
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

package nextflow.cli.pipeline

import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths

import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import groovy.transform.CompileStatic
import io.seqera.npr.api.schema.v1.CreatePipelineReleaseRequest
import io.seqera.npr.client.RegistryClient
import io.seqera.npr.client.RegistryException
import nextflow.cli.CmdBase
import nextflow.config.ConfigBuilder
import nextflow.config.RegistryConfig
import nextflow.exception.AbortOperationException
import nextflow.module.RegistryClientFactory
import nextflow.pipeline.GitProvenance
import nextflow.pipeline.PipelineBundle
import nextflow.pipeline.PipelineManager
import nextflow.util.TestOnly

/**
 * Pipeline publish subcommand
 */
@CompileStatic
@Parameters(commandDescription = "Publish a pipeline to the registry")
class CmdPipelinePublish extends CmdBase {

    static final long MAX_BUNDLE_SIZE = 50L * 1024 * 1024

    @Parameter(description = "[path]")
    List<String> args

    @Parameter(names = ["-version"], description = "Release version (default: manifest.version)")
    String version

    @Parameter(names = ["-force"], description = "Publish even if the Git checkout has uncommitted changes", arity = 0)
    boolean force = false

    @Parameter(names = ["-dry-run"], description = "Validate without uploading", arity = 0)
    boolean dryRun = false

    @Parameter(names = ["-registry"], description = "Target registry URL.")
    String registryUrl

    @TestOnly
    protected RegistryClient client

    @Override
    String getName() {
        return 'publish'
    }

    @Override
    void run() {
        if( args && args.size() > 1 )
            throw new AbortOperationException("Unexpected number of parameters")

        final dir = Paths.get(args ? args[0] : '.').toAbsolutePath().normalize()
        if( !Files.isDirectory(dir) )
            throw new AbortOperationException("Pipeline directory not found: ${dir}")

        final config = new ConfigBuilder().setBaseDir(dir).build(configFiles(dir))
        final manifest = (config.get('manifest') ?: [:]) as Map
        final name = manifest.name as String
        final ver = version ?: manifest.version as String
        if( !name || !PipelineManager.NAME_PATTERN.matcher(name).matches() )
            throw new AbortOperationException("Pipeline name must be set as 'manifest.name' in the form 'scope/name' (found: ${name ?: 'none'})")
        if( !ver || !(ver ==~ /\w[\w.+-]*/) )
            throw new AbortOperationException("Pipeline version must be set as 'manifest.version' or with -version (found: ${ver ?: 'none'})")

        if( !PipelineBundle.files(dir) )
            throw new AbortOperationException("No files would be published from ${dir} -- note that untracked files are not included in a Git checkout")

        final git = GitProvenance.of(dir)
        if( git && !git.clean && !force )
            throw new AbortOperationException("The Git checkout has uncommitted changes -- commit them or use -force to publish anyway")

        final bundle = Files.createTempFile("nf-pipeline-publish-", ".tar.gz")
        try {
            final checksum = PipelineBundle.create(dir, bundle)
            final size = Files.size(bundle)
            if( size > MAX_BUNDLE_SIZE )
                throw new AbortOperationException("The pipeline bundle is ${String.format('%.1f', size / (1024d * 1024d))} MB, over the 50 MB limit")

            final request = new CreatePipelineReleaseRequest()
                .version(ver)
                .checksum(checksum)
                .size(size)
                .description(manifest.description as String)
                .homePage(manifest.homePage as String)
            if( git?.repositoryUrl ) {
                request.repositoryUrl(git.repositoryUrl).gitRevision(git.revision).commitId(git.commitId)
            }
            else if( git ) {
                println "Note: the Git origin has no public URL, publishing without Git provenance"
            }

            if( dryRun ) {
                printDryRun(name, request)
                return
            }
            publish(config, name, request, bundle)
        }
        finally {
            Files.deleteIfExists(bundle)
        }
    }

    private void publish(Map config, String name, CreatePipelineReleaseRequest request, Path bundle) {
        final registry = client ?: newClient(config)
        println "Publishing pipeline ${name}@${request.version}..."
        try {
            final releaseId = registry.createPipelineRelease(name, request).releaseId
            try {
                registry.uploadPipelineRelease(releaseId, bundle)
            }
            catch( RegistryException e ) {
                throw new AbortOperationException("Upload of ${name}@${request.version} failed: ${e.message} -- run the command again to retry; the draft expires in 1 hour", e)
            }
        }
        catch( RegistryException e ) {
            throw new AbortOperationException("Cannot publish ${name}@${request.version}: ${e.message}", e)
        }
        println "Published pipeline ${name}@${request.version}"
    }

    private void printDryRun(String name, CreatePipelineReleaseRequest request) {
        println "Pipeline details:"
        println "  Name: ${name}"
        println "  Version: ${request.version}"
        println "  Size: ${request.size} bytes"
        println "  Checksum: ${request.checksum}"
        if( request.repositoryUrl )
            println "  Git: ${request.repositoryUrl} @ ${request.gitRevision} (${request.commitId})"
        println ""
        println "Dry run complete. Run without -dry-run to publish to the registry."
    }

    protected RegistryClient newClient(Map config) {
        return RegistryClientFactory.forConfig(registryConfig(config))
    }

    protected RegistryConfig registryConfig(Map config) {
        final registry = config.get('registry') as Map
        final opts = new LinkedHashMap(registry ?: [:])
        if( registryUrl )
            opts.url = registryUrl
        return new RegistryConfig(opts)
    }

    private static List<Path> configFiles(Path dir) {
        final file = dir.resolve('nextflow.config')
        return Files.exists(file) ? [file] : []
    }
}
