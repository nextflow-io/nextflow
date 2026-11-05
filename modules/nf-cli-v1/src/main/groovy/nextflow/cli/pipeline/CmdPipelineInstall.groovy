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

import java.nio.file.Path
import java.nio.file.Paths

import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import groovy.transform.CompileStatic
import io.seqera.npr.client.RegistryClient
import nextflow.cli.CmdBase
import nextflow.exception.AbortOperationException
import nextflow.module.RegistryClientFactory
import nextflow.pipeline.PipelineGitFetcher
import nextflow.pipeline.PipelineManager
import nextflow.util.TestOnly

/**
 * Pipeline install subcommand
 */
@CompileStatic
@Parameters(commandDescription = "Install a pipeline from the registry")
class CmdPipelineInstall extends CmdBase {

    @Parameter(description = "<pipeline>", required = true)
    List<String> args

    @Parameter(names = ["-version"], description = "Pipeline version (default: latest)")
    String version

    @Parameter(names = ["-force"], description = "Overwrite local modifications", arity = 0)
    boolean force = false

    @TestOnly
    protected RegistryClient client

    @TestOnly
    protected PipelineGitFetcher fetcher

    @TestOnly
    protected Path root

    @Override
    String getName() {
        return 'install'
    }

    @Override
    void run() {
        if( !args || args.size() != 1 )
            throw new AbortOperationException("Unexpected number of parameters")

        final baseDir = root ?: Paths.get('.').toAbsolutePath().normalize()
        final client = this.client ?: RegistryClientFactory.forConfig(PipelineCmdSupport.registryConfig(this, baseDir))
        final manager = fetcher ? new PipelineManager(baseDir, client, fetcher) : new PipelineManager(baseDir, client)
        final result = manager.install(args[0], version, force)
        if( result.upToDate() )
            println "Pipeline ${result.info().name}@${result.info().version} is already installed and up to date"
    }
}
