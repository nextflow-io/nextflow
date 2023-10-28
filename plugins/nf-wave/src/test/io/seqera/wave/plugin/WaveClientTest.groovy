/*
 * Copyright 2013-2023, Seqera Labs
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

package io.seqera.wave.plugin

import static java.nio.file.StandardOpenOption.*

import java.net.http.HttpRequest
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.attribute.FileTime

import com.sun.net.httpserver.HttpExchange
import com.sun.net.httpserver.HttpHandler
import com.sun.net.httpserver.HttpServer
import groovy.json.JsonOutput
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.SysEnv
import nextflow.container.inspect.ContainerInspectMode
import nextflow.container.inspect.ContainersInspector
import nextflow.extension.FilesEx
import nextflow.file.FileHelper
import nextflow.processor.TaskRun
import nextflow.script.bundle.ResourcesBundle
import org.apache.commons.compress.archivers.ArchiveStreamFactory
import org.apache.commons.compress.archivers.tar.TarArchiveEntry
import org.apache.commons.compress.archivers.tar.TarArchiveInputStream
import org.apache.commons.compress.compressors.gzip.GzipCompressorInputStream
import spock.lang.Shared
import spock.lang.Specification
import spock.lang.Unroll
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class WaveClientTest extends Specification {

    @CompileStatic
    private static List<Path> untar(final InputStream is, final Path outputDir)  {

        final untaredFiles = new LinkedList<Path>();
        final debInputStream = (TarArchiveInputStream) new ArchiveStreamFactory().createArchiveInputStream("tar", is);
        TarArchiveEntry entry = null;
        while ((entry = (TarArchiveEntry)debInputStream.getNextEntry()) != null) {
            final outputFile = outputDir.resolve(entry.getName())
            if (entry.isDirectory()) {
                if (!outputFile.exists()) {
                    if (!outputFile.mkdirs()) {
                        throw new IllegalStateException(String.format("Couldn't create directory %s.", outputFile));
                    }
                    outputFile.setPermissionsMode( entry.getMode() )
                    outputFile.setLastModified( entry.getLastModifiedDate().getTime() )
                }
            }
            else {
                println "outputFile=$outputFile; mode=${entry.mode}"
                outputFile.parent.mkdirs()
                final outputFileStream = Files.newOutputStream(outputFile, CREATE, APPEND)
                debInputStream.transferTo(outputFileStream)
                outputFileStream.close()
                outputFile.setPermissionsMode(entry.getMode())
                outputFile.setLastModified( entry.getLastModifiedDate().getTime() )
            }
            untaredFiles.add(outputFile);
        }
        debInputStream.close();

        return untaredFiles
    }

    byte[] uncompress( byte[] bytes ) {
        try (def stream = new GzipCompressorInputStream(new ByteArrayInputStream(bytes))) {
            def buffer = new ByteArrayOutputStream()
            stream.transferTo(buffer)
            return buffer.toByteArray()
        }
    }

    def 'should tar file' () {
        given:
        def LAST_MODIFIED = FileTime.fromMillis(1_000_000_000_000)
        def sess = Mock(Session) { getConfig() >> [wave:[:]]}
        def folder = Files.createTempDirectory('test')
        and:
        def bundlePath = folder.resolve('bundle'); bundlePath.mkdir()
        bundlePath.resolve('main.nf').text = "I'm the main file"
        bundlePath.resolve('this/that').mkdirs()
        Files.write(bundlePath.resolve('this/hola.txt'), "Hola".bytes)
        Files.write(bundlePath.resolve('this/hello.txt'), "Hello".bytes)
        Files.write(bundlePath.resolve('this/that/ciao.txt'), "Ciao".bytes)
        and:
        FileHelper.visitFiles([type:'any'], bundlePath, '**', {
            Files.setLastModifiedTime(it, LAST_MODIFIED)
            final mode = it.isDirectory() ? 0700 : 0600
            FilesEx.setPermissionsMode(it, mode)
        })
        and:
        def wave = new WaveClient(sess)

        when:
        def bundle = ResourcesBundle.scan(bundlePath)
        def layer = wave.makeLayer(bundle)
        then:
        layer.tarDigest == 'sha256:81200f6ad32793567d8070375dc51312a1711fedf6a1c6f5e4a97fa3014f3491'
        layer.gzipDigest == 'sha256:09a2deca4293245909223db505cf69affa1a8ff8acb745fe3cad38bc0b719110'
        layer.gzipSize == 254
        and:
        def gzip = layer.location.replace('data:','').decodeBase64()
        def tar = uncompress(gzip)
        def result = folder.resolve('result')
        untar( new ByteArrayInputStream(tar), result)
        and:
        result.resolve('main.nf').text == bundlePath.resolve('main.nf').text
        result.resolve('this/hola.txt').text == bundlePath.resolve('this/hola.txt').text
        result.resolve('this/hello.txt').text == bundlePath.resolve('this/hello.txt').text
        result.resolve('this/that/ciao.txt').text == bundlePath.resolve('this/that/ciao.txt').text

        /*
         * create a bundle using different base directory
         */
        when:
        bundle = ResourcesBundle.scan(bundlePath, [baseDirectory: 'usr/local'])
        layer = wave.makeLayer(bundle)
        then:
        def gzip2 = layer.location.replace('data:','').decodeBase64()
        def tar2 = uncompress(gzip2)
        def result2 = folder.resolve('result2')
        untar( new ByteArrayInputStream(tar2), result2)
        and:
        result2.resolve('usr/local/main.nf').text == bundlePath.resolve('main.nf').text
        result2.resolve('usr/local/this/hola.txt').text == bundlePath.resolve('this/hola.txt').text
        result2.resolve('usr/local/this/hello.txt').text == bundlePath.resolve('this/hello.txt').text
        result2.resolve('usr/local/this/that/ciao.txt').text == bundlePath.resolve('this/that/ciao.txt').text

        cleanup:
        folder?.deleteDir()
    }

    def 'should create request object' () {
        given:
        def session = Mock(Session) { getConfig() >> [:]}
        def IMAGE =  'foo:latest'
        def wave = new WaveClient(session)

        when:
        def req = wave.makeRequest(WaveAssets.fromImage(IMAGE))
        then:
        req.containerImage == IMAGE
        !req.containerPlatform
        !req.containerFile
        !req.condaFile
        !req.spackFile
        !req.containerConfig.layers
        !req.freeze
        !req.dryRun
        and:
        req.fingerprint == 'bd2cb4b32df41f2d290ce2366609f2ad'
        req.timestamp instanceof String
    }

    def 'should create request object with freeze mode' () {
        given:
        def session = Mock(Session) { getConfig() >> [wave:[freeze:true]]}
        def IMAGE =  'foo:latest'
        def wave = new WaveClient(session)

        when:
        def req = wave.makeRequest(WaveAssets.fromImage(IMAGE))
        then:
        req.containerImage == IMAGE
        !req.containerPlatform
        !req.containerFile
        !req.condaFile
        !req.spackFile
        !req.containerConfig.layers
        and:
        req.freeze
        and:
        req.fingerprint == 'bd2cb4b32df41f2d290ce2366609f2ad'
        req.timestamp instanceof String
    }

    def 'should create request object with dry-run mode' () {
        given:
        ContainerInspectMode.activate(true)
        def session = Mock(Session) { getConfig() >> [:]}
        def IMAGE =  'foo:latest'
        def wave = new WaveClient(session)

        when:
        def req = wave.makeRequest(WaveAssets.fromImage(IMAGE))
        then:
        req.containerImage == IMAGE
        !req.containerPlatform
        !req.containerFile
        !req.condaFile
        !req.spackFile
        !req.containerConfig.layers
        and:
        req.dryRun
        and:
        req.fingerprint == 'bd2cb4b32df41f2d290ce2366609f2ad'
        req.timestamp instanceof String

        cleanup:
        ContainerInspectMode.activate(false)
    }

    def 'should create request object and platform' () {
        given:
        def session = Mock(Session) { getConfig() >> [:]}
        def IMAGE =  'foo:latest'
        def PLATFORM = 'amd64'
        def wave = new WaveClient(session)

        when:
        def req = wave.makeRequest(WaveAssets.fromImage(IMAGE, PLATFORM))
        then:
        req.containerImage == IMAGE
        req.containerPlatform == PLATFORM
        and:
        !req.containerFile
        !req.condaFile
        !req.spackFile
        !req.containerConfig.layers
        and:
        req.fingerprint == 'd31044e6594126479585c0cdca15c15e'
        req.timestamp instanceof String
    }

    def 'should create request object with dockerfile' () {
        given:
        def session = Mock(Session) { getConfig() >> [:]}
        def DOCKERFILE =  'FROM foo:latest\nRUN something'
        def wave = new WaveClient(session)

        when:
        def req = wave.makeRequest(WaveAssets.fromDockerfile(DOCKERFILE))
        then:
        !req.containerImage
        new String(req.containerFile.decodeBase64()) == DOCKERFILE
        !req.condaFile
        !req.spackFile
        !req.containerConfig.layers
    }

    def 'should create request object with singularityfile' () {
        given:
        def session = Mock(Session) { getConfig() >> [:]}
        def SINGULARITY_FILE =  'From foo:latest'
        def wave = new WaveClient(session)
        and:
        def assets = new WaveAssets(null,
                'linux/amd64',
                null,
                null,
                SINGULARITY_FILE,
                null,
                null,
                null,
                true)
        when:
        def req = wave.makeRequest(assets)
        then:
        !req.containerImage
        new String(req.containerFile.decodeBase64()) == SINGULARITY_FILE
        !req.condaFile
        !req.spackFile
        !req.containerConfig.layers
        and:
        req.format == 'sif'
    }

    def 'should create request object with build and cache repos' () {
        given:
        def session = Mock(Session) { getConfig() >> [wave:[build:[repository:'some/repo',cacheRepository:'some/cache']]]}
        def DOCKERFILE =  'FROM foo:latest\nRUN something'
        def wave = new WaveClient(session)

        when:
        def req = wave.makeRequest(WaveAssets.fromDockerfile(DOCKERFILE))
        then:
        req.buildRepository == 'some/repo'
        req.cacheRepository == 'some/cache'
        !req.containerImage
        new String(req.containerFile.decodeBase64()) == DOCKERFILE
        !req.condaFile
        !req.spackFile
        !req.containerConfig.layers
    }

    def 'should create request object with conda file' () {
        given:
        def folder = Files.createTempDirectory('test')
        and:
        def DOCKERFILE = 'from foo:latest'
        def CONDAFILE = folder.resolve('conda.yml'); CONDAFILE.text = 'some conda recipe here'
        and:
        def session = Mock(Session) { getConfig() >> [:]}
        def wave = new WaveClient(session)

        when:
        def req = wave.makeRequest(new WaveAssets(null, null, null, null, DOCKERFILE, CONDAFILE))
        then:
        !req.containerImage
        new String(req.containerFile.decodeBase64()) == DOCKERFILE
        new String(req.condaFile.decodeBase64()) == CONDAFILE.text
        !req.containerConfig.layers

        cleanup:
        folder?.deleteDir()
    }

    def 'should create request object with spack file' () {
        given:
        def folder = Files.createTempDirectory('test')
        and:
        def DOCKERFILE = 'from foo:latest'
        def SPACKFILE = folder.resolve('spack.yaml'); SPACKFILE.text = 'some spack recipe here'
        and:
        def session = Mock(Session) { getConfig() >> [:]}
        def wave = new WaveClient(session)

        when:
        def req = wave.makeRequest(new WaveAssets(null, null, null, null, DOCKERFILE, null, SPACKFILE))
        then:
        !req.containerImage
        new String(req.containerFile.decodeBase64()) == DOCKERFILE
        new String(req.spackFile.decodeBase64()) == SPACKFILE.text
        !req.containerConfig.layers

        cleanup:
        folder?.deleteDir()
    }

    def 'should  create request with module resources' () {
        given:
        def MODULE_RES = Mock(ResourcesBundle) {hasEntries() >> true }
        def MODULE_LAYER = Mock(ContainerLayer)

        and:
        def session = Mock(Session) { getConfig() >> [:]}
        WaveClient wave = Spy(WaveClient, constructorArgs: [session])

        when:
        def assets = new WaveAssets('my:image', null, MODULE_RES)
        def req = wave.makeRequest(assets)
        then:
        1 * wave.makeLayer(MODULE_RES) >> MODULE_LAYER
        and:
        req.containerImage == 'my:image'
        req.containerConfig.layers.size()==1
        req.containerConfig.layers[0] == MODULE_LAYER
    }

    def 'should  create request with module and project resources' () {
        given:
        def MODULE_RES = Mock(ResourcesBundle) {hasEntries() >> true }
        def MODULE_LAYER = Mock(ContainerLayer)
        and:
        def PROJECT_RES = Mock(ResourcesBundle) { hasEntries() >> true }
        def PROJECT_LAYER = Mock(ContainerLayer)
        and:
        def session = Mock(Session) { getConfig() >> [:]}
        WaveClient wave = Spy(WaveClient, constructorArgs: [session])

        when:
        def assets = new WaveAssets('my:image', null, MODULE_RES, null, null, null, null, PROJECT_RES)
        def req = wave.makeRequest(assets)
        then:
        1 * wave.makeLayer(MODULE_RES) >> MODULE_LAYER
        1 * wave.makeLayer(PROJECT_RES) >> PROJECT_LAYER
        and:
        req.containerImage == 'my:image'
        req.containerConfig.layers.size()==2
        req.containerConfig.layers[0] == PROJECT_LAYER
        req.containerConfig.layers[1] == MODULE_LAYER
    }


    def 'should create asset with image' () {
        given:
        def session = Mock(Session) { getConfig() >> [:]}
        def task = Mock(TaskRun) { getConfig() >> [arch:'amd64'] }
        def IMAGE = 'foo:latest'
        and:
        def client = new WaveClient(session)

        when:
        def assets = client.resolveAssets(task, IMAGE, false)
        then:
        assets.containerImage == IMAGE
        !assets.moduleResources
        !assets.containerFile
        !assets.containerConfig
        !assets.condaFile
        !assets.spackFile
        !assets.projectResources
        assets.containerPlatform == 'linux/amd64'
    }

    def 'should create asset with image and platform' () {
        given:
        def ARCH = 'linux/arm64'
        def session = Mock(Session) { getConfig() >> [:] }
        def task = Mock(TaskRun) { getConfig() >> [arch:ARCH.toString()] }
        def IMAGE = 'foo:latest'
        and:
        def client = new WaveClient(session)

        when:
        def assets = client.resolveAssets(task, IMAGE, false)
        then:
        assets.containerImage == IMAGE
        assets.containerPlatform == 'linux/arm64'
        !assets.moduleResources
        !assets.containerFile
        !assets.containerConfig
        !assets.condaFile
        !assets.spackFile
        !assets.projectResources
    }

    def 'should create asset with image and bundle' () {
        given:
        def IMAGE = 'foo:latest'
        def BUNDLE = Mock(ResourcesBundle)
        and:
        def session = Mock(Session) { getConfig() >> [:]}
        def task = Mock(TaskRun) { getConfig() >> [:]; getModuleBundle() >> BUNDLE }
        and:
        def client = new WaveClient(session)

        when:
        def assets = client.resolveAssets(task, IMAGE, false)
        then:
        assets.containerImage == IMAGE
        assets.moduleResources == BUNDLE
        !assets.containerFile
        !assets.containerConfig
        !assets.condaFile
        !assets.spackFile
        !assets.projectResources
    }

    def 'should create asset with image and bundle and container config' () {
        given:
        def IMAGE = 'foo:latest'
        def BUNDLE = Mock(ResourcesBundle)
        def ARCH = 'linux/arm64'
        def CONTAINER_CONFIG = new ContainerConfig(entrypoint: ['entry.sh'], layers: [new ContainerLayer(location: 'http://somewhere')])
        and:
        def session = Mock(Session) { getConfig() >> [:]}
        def task = Mock(TaskRun) { getConfig() >> [arch:ARCH.toString()]; getModuleBundle() >> BUNDLE }
        and:
        WaveClient client = Spy(WaveClient, constructorArgs:[session])

        when:
        def assets = client.resolveAssets(task, IMAGE, false)
        then:
        client.resolveContainerConfig(ARCH) >> CONTAINER_CONFIG
        and:
        assets.containerImage == IMAGE
        assets.moduleResources == BUNDLE
        assets.containerConfig == CONTAINER_CONFIG
        and:
        !assets.containerFile
        !assets.condaFile
        !assets.spackFile
        !assets.projectResources
    }

    def 'should create asset with dockerfile' () {
        given:
        def folder = Files.createTempDirectory('test')
        def session = Mock(Session) { getConfig() >> [:]}
        and:
        def DOCKERFILE = folder.resolve('Dockerfile')
        DOCKERFILE.text = 'FROM foo\nRUN this/that'
        and:
        def BUNDLE = Mock(ResourcesBundle) { getDockerfile() >> DOCKERFILE }
        and:
        def task = Mock(TaskRun) {getModuleBundle() >> BUNDLE; getConfig() >> [:] }
        and:
        def client = new WaveClient(session)

        when:
        def assets = client.resolveAssets(task, null, false)
        then:
        assets.containerFile == 'FROM foo\nRUN this/that'
        assets.moduleResources == BUNDLE
        !assets.containerImage
        !assets.containerConfig
        !assets.condaFile
        !assets.spackFile
        !assets.projectResources

        cleanup:
        folder?.deleteDir()
    }

    def 'should create asset with conda recipe' () {
        given:
        def session = Mock(Session) { getConfig() >> [:]}
        and:
        def task = Mock(TaskRun) {getConfig() >> [conda:"bioconda::rseqc=3.0.1 'conda-forge::r-base>=3.5'"] }
        and:
        def client = new WaveClient(session)

        when:
        def assets = client.resolveAssets(task, null, false)
        then:
        assets.containerFile == '''\
                FROM mambaorg/micromamba:1.5.1
                COPY --chown=$MAMBA_USER:$MAMBA_USER conda.yml /tmp/conda.yml
                RUN micromamba install -y -n base -f /tmp/conda.yml \\
                    && micromamba install -y -n base conda-forge::procps-ng \\
                    && micromamba clean -a -y
                USER root
                    '''.stripIndent()
        and:
        !assets.moduleResources
        !assets.containerImage
        !assets.containerConfig
        !assets.spackFile
        !assets.projectResources
        and:
        assets.condaFile.text == '''\
                channels:
                - seqera
                - conda-forge
                - bioconda
                - defaults
                dependencies:
                - bioconda::rseqc=3.0.1
                - conda-forge::r-base>=3.5
                '''.stripIndent(true)
    }

    def 'should create asset with conda lock file' () {
        given:
        def session = Mock(Session) { getConfig() >> [:]}
        and:
        def task = Mock(TaskRun) {getConfig() >> [conda:'https://host.com/conda-lock.yml'] }
        and:
        def client = new WaveClient(session)

        when:
        def assets = client.resolveAssets(task, null, false)
        then:
        assets.containerFile == '''\
                FROM mambaorg/micromamba:1.5.1
                RUN \\
                    micromamba install -y -n base -c seqera -c conda-forge -c bioconda -c defaults -f https://host.com/conda-lock.yml \\
                    && micromamba install -y -n base conda-forge::procps-ng \\
                    && micromamba clean -a -y
                USER root
                    '''.stripIndent()
        and:
        !assets.moduleResources
        !assets.containerImage
        !assets.containerConfig
        !assets.condaFile
        !assets.spackFile
        !assets.projectResources
    }

    def 'should create asset with spack recipe' () {
        given:
        def session = Mock(Session) { getConfig() >> [:]}
        and:
        def task = Mock(TaskRun) {getConfig() >> [spack:"rseqc@3.0.1 'rbase@3.5'", arch:"amd64"] }
        and:
        def client = new WaveClient(session)

        when:
        def assets = client.resolveAssets(task, null, false)
        then:
        assets.containerFile == '''\
                # Runner image
                FROM {{spack_runner_image}}
                
                COPY --from=builder /opt/spack-env /opt/spack-env
                COPY --from=builder /opt/software /opt/software
                COPY --from=builder /opt/._view /opt/._view
                
                # Entrypoint for Singularity
                RUN mkdir -p /.singularity.d/env && \\
                    cp -p /opt/spack-env/z10_spack_environment.sh /.singularity.d/env/91-environment.sh
                # Entrypoint for Docker
                RUN echo "#!/usr/bin/env bash\\n\\nset -ef -o pipefail\\nsource /opt/spack-env/z10_spack_environment.sh\\nexec \\"\\\$@\\"" \\
                    >/opt/spack-env/spack_docker_entrypoint.sh && chmod a+x /opt/spack-env/spack_docker_entrypoint.sh
                
                
                ENTRYPOINT [ "/opt/spack-env/spack_docker_entrypoint.sh" ]
                CMD [ "/bin/bash" ]
                '''.stripIndent()

        and:
        !assets.moduleResources
        !assets.containerImage
        !assets.containerConfig
        !assets.condaFile
        !assets.projectResources
        and:
        assets.spackFile.text == '''\
                spack:
                  specs: [rseqc@3.0.1, rbase@3.5]
                  concretizer: {unify: true, reuse: false}
                '''.stripIndent(true)
    }

    def 'should create asset with conda file' () {
        given:
        def folder = Files.createTempDirectory('test')
        def condaFile = folder.resolve('conda.yml'); condaFile.text = 'the-conda-recipe-here'
        and:
        def session = Mock(Session) { getConfig() >> [:]}
        def task = Mock(TaskRun) {getConfig() >> [conda:condaFile.toString()] }
        and:
        def client = new WaveClient(session)

        when:
        def assets = client.resolveAssets(task, null, false)
        then:
        assets.containerFile == '''\
                FROM mambaorg/micromamba:1.5.1
                COPY --chown=$MAMBA_USER:$MAMBA_USER conda.yml /tmp/conda.yml
                RUN micromamba install -y -n base -f /tmp/conda.yml \\
                    && micromamba install -y -n base conda-forge::procps-ng \\
                    && micromamba clean -a -y
                USER root
                    '''.stripIndent()
        and:
        assets.condaFile == condaFile
        and:
        !assets.moduleResources
        !assets.containerImage
        !assets.containerConfig
        !assets.spackFile
        !assets.projectResources

        cleanup:
        folder?.deleteDir()
    }

    def 'should create asset with spack file' () {
        given:
        def folder = Files.createTempDirectory('test')
        def spackFile = folder.resolve('spack.yaml'); spackFile.text = 'the-spack-recipe-here'
        and:
        def session = Mock(Session) { getConfig() >> [:]}
        def task = Mock(TaskRun) {getConfig() >> [spack:spackFile.toString(), arch: 'amd64'] }
        and:
        def client = new WaveClient(session)

        when:
        def assets = client.resolveAssets(task, null, false)
        then:
        assets.containerFile == '''\
                # Runner image
                FROM {{spack_runner_image}}
                
                COPY --from=builder /opt/spack-env /opt/spack-env
                COPY --from=builder /opt/software /opt/software
                COPY --from=builder /opt/._view /opt/._view
                
                # Entrypoint for Singularity
                RUN mkdir -p /.singularity.d/env && \\
                    cp -p /opt/spack-env/z10_spack_environment.sh /.singularity.d/env/91-environment.sh
                # Entrypoint for Docker
                RUN echo "#!/usr/bin/env bash\\n\\nset -ef -o pipefail\\nsource /opt/spack-env/z10_spack_environment.sh\\nexec \\"\\\$@\\"" \\
                    >/opt/spack-env/spack_docker_entrypoint.sh && chmod a+x /opt/spack-env/spack_docker_entrypoint.sh
                
                
                ENTRYPOINT [ "/opt/spack-env/spack_docker_entrypoint.sh" ]
                CMD [ "/bin/bash" ]
                '''.stripIndent()
        and:
        assets.spackFile == spackFile
        and:
        !assets.moduleResources
        !assets.containerImage
        !assets.containerConfig
        !assets.condaFile
        !assets.projectResources

        cleanup:
        folder?.deleteDir()
    }

    // ==== singularity native build + conda ====

    def 'should create asset with conda recipe and singularity native build' () {
        given:
        def session = Mock(Session) { getConfig() >> [:]}
        and:
        def task = Mock(TaskRun) {getConfig() >> [conda:'salmon=1.2.3'] }
        and:
        def client = new WaveClient(session)

        when:
        def assets = client.resolveAssets(task, null, true)
        then:
        assets.containerFile == '''\
                BootStrap: docker
                From: mambaorg/micromamba:1.5.1
                %files
                    {{wave_context_dir}}/conda.yml /scratch/conda.yml
                %post
                    micromamba install -y -n base -f /scratch/conda.yml
                    micromamba install -y -n base conda-forge::procps-ng
                    micromamba clean -a -y
                %environment
                    export PATH="$MAMBA_ROOT_PREFIX/bin:$PATH"
                    '''.stripIndent()
        and:
        assets.singularity
        and:
        !assets.moduleResources
        !assets.containerImage
        !assets.containerConfig
        !assets.spackFile
        !assets.projectResources
        and:
        assets.condaFile.text == '''\
                channels:
                - seqera
                - conda-forge
                - bioconda
                - defaults
                dependencies:
                - salmon=1.2.3
                '''.stripIndent(true)
    }

    def 'should create asset with conda remote lock file and singularity native build' () {
        given:
        def session = Mock(Session) { getConfig() >> [:]}
        and:
        def task = Mock(TaskRun) {getConfig() >> [conda:'https://host.com/lock-file.yaml'] }
        and:
        def client = new WaveClient(session)

        when:
        def assets = client.resolveAssets(task, null, true)
        then:
        assets.containerFile == '''\
                BootStrap: docker
                From: mambaorg/micromamba:1.5.1
                %post
                    micromamba install -y -n base -c seqera -c conda-forge -c bioconda -c defaults -f https://host.com/lock-file.yaml
                    micromamba install -y -n base conda-forge::procps-ng
                    micromamba clean -a -y
                %environment
                    export PATH="$MAMBA_ROOT_PREFIX/bin:$PATH"
                    '''.stripIndent()
        and:
        assets.singularity
        and:
        !assets.moduleResources
        !assets.containerImage
        !assets.containerConfig
        !assets.condaFile
        !assets.spackFile
        !assets.projectResources
    }

    def 'should create asset with conda file and singularity native build' () {
        given:
        def folder = Files.createTempDirectory('test')
        def condaFile = folder.resolve('conda.yml'); condaFile.text = 'the-conda-recipe-here'
        and:
        def session = Mock(Session) { getConfig() >> [:]}
        def task = Mock(TaskRun) {getConfig() >> [conda:condaFile.toString()] }
        and:
        def client = new WaveClient(session)

        when:
        def assets = client.resolveAssets(task, null, true)
        then:
        assets.containerFile == '''\
                BootStrap: docker
                From: mambaorg/micromamba:1.5.1
                %files
                    {{wave_context_dir}}/conda.yml /scratch/conda.yml
                %post
                    micromamba install -y -n base -f /scratch/conda.yml
                    micromamba install -y -n base conda-forge::procps-ng
                    micromamba clean -a -y
                %environment
                    export PATH="$MAMBA_ROOT_PREFIX/bin:$PATH"                    
                '''.stripIndent()
        and:
        assets.condaFile == condaFile
        assets.singularity
        and:
        !assets.moduleResources
        !assets.containerImage
        !assets.containerConfig
        !assets.spackFile
        !assets.projectResources

        cleanup:
        folder?.deleteDir()
    }

    def 'should create assets with spack recipe for singularity' () {
        given:
        def session = Mock(Session) { getConfig() >> [wave:[build:[spack:[commands: ['cmd-foo','cmd-bar']]]]]}
        and:
        def task = Mock(TaskRun) {getConfig() >> [spack:"rseqc@3.0.1 'rbase@3.5'", arch:"amd64"] }
        and:
        def client = new WaveClient(session)

        when:
        def assets = client.resolveAssets(task, null, true)
        then:
        assets.containerFile == '''\
                Bootstrap: docker
                From: {{spack_runner_image}}
                stage: final
                 
                %files from build
                    /opt/spack-env /opt/spack-env
                    /opt/software /opt/software
                    /opt/._view /opt/._view
                    /opt/spack-env/z10_spack_environment.sh /.singularity.d/env/91-environment.sh
                 
                %post
                    cmd-foo
                    cmd-bar
                '''.stripIndent()

        and:
        !assets.moduleResources
        !assets.containerImage
        !assets.containerConfig
        !assets.condaFile
        !assets.projectResources
        and:
        assets.spackFile.text == '''\
                spack:
                  specs: [rseqc@3.0.1, rbase@3.5]
                  concretizer: {unify: true, reuse: false}
                '''.stripIndent(true)
    }

    def 'should create asset with spack file for singularity' () {
        given:
        def folder = Files.createTempDirectory('test')
        def spackFile = folder.resolve('spack.yml');
        spackFile.text = '''\
                spack:
                  specs: [rseqc@3.0.1, rbase@3.5]
                  concretizer: {unify: true, reuse: false}
                '''.stripIndent(true)
        and:
        def session = Mock(Session) { getConfig() >> [wave:[build:[spack:[basePackages: 'nano@1.2.3']]]]}
        def task = Mock(TaskRun) {getConfig() >> [spack:spackFile.toString()] }
        and:
        def client = new WaveClient(session)

        when:
        def assets = client.resolveAssets(task, null, true)
        then:
        assets.containerFile == '''\
                    Bootstrap: docker
                    From: {{spack_runner_image}}
                    stage: final
                     
                    %files from build
                        /opt/spack-env /opt/spack-env
                        /opt/software /opt/software
                        /opt/._view /opt/._view
                        /opt/spack-env/z10_spack_environment.sh /.singularity.d/env/91-environment.sh
                     
                    %post
                    '''.stripIndent()
        and:
        !assets.moduleResources
        !assets.containerImage
        !assets.containerConfig
        !assets.projectResources
        !assets.condaFile
        and:
        assets.spackFile.text == '''\
                spack:
                  specs: [rseqc@3.0.1, rbase@3.5, nano@1.2.3]
                  concretizer: {unify: true, reuse: false}
                '''.stripIndent(true)

        cleanup:
        folder?.deleteDir()
    }

    def 'should create assets with project resources' () {
        given:
        def MODULE_RES = Mock(ResourcesBundle)
        def PROJECT_RES = Mock(ResourcesBundle)
        def CONTAINER_CONFIG = Mock(ContainerConfig)
        def BIN_DIR = Path.of('/something/bin')
        def ARCH = 'linux/arm64'
        and:
        def task = Mock(TaskRun) {getModuleBundle() >> MODULE_RES; getConfig() >> [arch:ARCH.toString()] }
        and:
        def session = Mock(Session) {
            getConfig() >> [wave: [bundleProjectResources: true]]
            getBinDir() >> BIN_DIR
        }
        and:
        WaveClient wave = Spy(WaveClient, constructorArgs: [session])

        when:
        def assets = wave.resolveAssets(task, 'image:latest', false)
        then:
        1 * wave.projectResources(BIN_DIR) >> PROJECT_RES
        and:
        1 * wave.resolveContainerConfig(ARCH) >> CONTAINER_CONFIG
        and:
        assets.moduleResources == MODULE_RES
        assets.projectResources ==  PROJECT_RES
        assets.containerImage == 'image:latest'
    }

    def 'should resolve conflicts' () {
        given:
        def session = Mock(Session) { getConfig() >> [:]}
        and:
        def client = new WaveClient(session)

        when:
        def result = client.resolveConflicts([:], [])
        then:
        result == [:]

        when:
        result = client.resolveConflicts([dockerfile:'x',conda:'y',container:'z'], ['conda'])
        then:
        result == [conda:'y']

        when:
        result = client.resolveConflicts([dockerfile:'x',conda:'y',container:'z'], ['conda','dockerfile'])
        then:
        result == [conda:'y']

        when:
        result = client.resolveConflicts([dockerfile:'x',container:'z'], ['conda','dockerfile'])
        then:
        result == [dockerfile:'x']

        when:
        result = client.resolveConflicts([spack:'x',container:'z'], ['conda','spack'])
        then:
        result == [spack:'x']
    }

    def 'should patch strategy for singularity' () {
        given:
        def session = Mock(Session) { getConfig() >> [:]}
        and:
        def client = new WaveClient(session)

        expect:
        client.patchStrategy(Collections.unmodifiableList(STRATEGY), SING) == EXPECTED

        where:
        STRATEGY                                            | SING      | EXPECTED
        ['conda','dockerfile', 'spack']                     | false     | ['conda','dockerfile', 'spack']
        ['conda','dockerfile', 'spack']                     | true      | ['conda','singularityfile', 'spack']
        ['conda','dockerfile', 'spack']                     | true      | ['conda','singularityfile', 'spack']
        ['conda','singularityfile','dockerfile', 'spack']   | true      | ['conda','singularityfile','dockerfile', 'spack']
    }

    def 'should check conflicts' () {
        given:
        def session = Mock(Session) { getConfig() >> [:]}
        and:
        def client = new WaveClient(session)

        when:
        client.checkConflicts([:], 'foo')
        then:
        noExceptionThrown()

        when:
        client.checkConflicts([conda:'this', container:'that'], 'foo')
        then:
        def e = thrown(IllegalArgumentException)
        e.message == "Process 'foo' declares both 'container' and 'conda' directives that conflict each other"

        when:
        client.checkConflicts([conda:'this', dockerfile:'that'], 'foo')
        then:
        e = thrown(IllegalArgumentException)
        e.message == "Process 'foo' declares both a 'conda' directive and a module bundle dockerfile that conflict each other"

        when:
        client.checkConflicts([container:'this', dockerfile:'that'], 'foo')
        then:
        e = thrown(IllegalArgumentException)
        e.message == "Process 'foo' declares both a 'container' directive and a module bundle dockerfile that conflict each other"

        when:
        client.checkConflicts([spack:'this', dockerfile:'that'], 'foo')
        then:
        e = thrown(IllegalArgumentException)
        e.message == "Process 'foo' declares both a 'spack' directive and a module bundle dockerfile that conflict each other"

        when:
        client.checkConflicts([spack:'this', container:'that'], 'foo')
        then:
        e = thrown(IllegalArgumentException)
        e.message == "Process 'foo' declares both 'container' and 'spack' directives that conflict each other"

        when:
        client.checkConflicts([conda:'this', spack:'that'], 'foo')
        then:
        e = thrown(IllegalArgumentException)
        e.message == "Process 'foo' declares both 'spack' and 'conda' directives that conflict each other"

        // singularity file checks
        when:
        client.checkConflicts([conda:'this', singularityfile:'that'], 'foo')
        then:
        e = thrown(IllegalArgumentException)
        e.message == "Process 'foo' declares both a 'conda' directive and a module bundle singularityfile that conflict each other"

        when:
        client.checkConflicts([container:'this', singularityfile:'that'], 'foo')
        then:
        e = thrown(IllegalArgumentException)
        e.message == "Process 'foo' declares both a 'container' directive and a module bundle singularityfile that conflict each other"

        when:
        client.checkConflicts([spack:'this', singularityfile:'that'], 'foo')
        then:
        e = thrown(IllegalArgumentException)
        e.message == "Process 'foo' declares both a 'spack' directive and a module bundle singularityfile that conflict each other"

    }

    def 'should get project resource bundle' () {
        given:
        def folder = Files.createTempDirectory('test')
        def root = folder.resolve('bundle'); root.mkdir()
        root.resolve('main.nf').text = "I'm the main file"
        root.resolve('bin/nested').mkdirs()
        root.resolve('this/that').mkdirs()
        Files.write(root.resolve('bin/hola.sh'), "Hola".bytes)
        Files.write(root.resolve('bin/hello.sh'), "Hello".bytes)
        Files.write(root.resolve('bin/nested/script.sh'), "Script".bytes)
        Files.write(root.resolve('this/that/ciao.txt'), "Ciao".bytes)
        and:
        def sess = Mock(Session) {getConfig() >> [wave:[:]] }
        and:
        def wave = new WaveClient(sess)

        when:
        def result = wave.projectResources(root.resolve('bin'))
        then:
        result.getEntries() ==  ['usr/local/bin/hola.sh','usr/local/bin/hello.sh','usr/local/bin/nested','usr/local/bin/nested/script.sh'] as Set

        expect:
        wave.projectResources(null) == null

        cleanup:
        folder?.deleteDir()

    }

    def 'should send request with tower access token' () {
        given:
        def config = [wave:[:], tower:[accessToken:'foo', workspaceId:123, endpoint: 'http://foo.com']]
        def sess = Mock(Session) {getConfig() >> config }
        and:
        def wave = Spy(new WaveClient(sess))
        def assets = Mock(WaveAssets)
        def request = new SubmitContainerTokenRequest()
        when:
        wave.sendRequest(assets)
        then:
        1 * wave.makeRequest(assets) >> request
        1 * wave.sendRequest(request) >> { List it ->
            assert (it[0] == request)
            assert (it[0] as SubmitContainerTokenRequest).towerAccessToken == 'foo'
            assert (it[0] as SubmitContainerTokenRequest).towerWorkspaceId == 123
            assert (it[0] as SubmitContainerTokenRequest).towerEndpoint == 'http://foo.com'
        }
    }

    def 'should send request with tower access token and refresh token' () {
        given:
        SysEnv.push([TOWER_WORKFLOW_ID:'1234', TOWER_REFRESH_TOKEN: 'xyz', TOWER_ACCESS_TOKEN: 'foo'])
        and:
        def config = [wave:[:], tower:[endpoint: 'http://foo.com']]
        def sess = Mock(Session) {getConfig() >> config }
        and:
        def wave = Spy(new WaveClient(sess))
        def assets = Mock(WaveAssets)
        def request = new SubmitContainerTokenRequest()
        when:
        wave.sendRequest(assets)
        then:
        1 * wave.makeRequest(assets) >> request
        1 * wave.sendRequest(request) >> { List it ->
            assert (it[0] == request)
            assert (it[0] as SubmitContainerTokenRequest).towerAccessToken == 'foo'
            assert (it[0] as SubmitContainerTokenRequest).towerRefreshToken == 'xyz'
            assert (it[0] as SubmitContainerTokenRequest).towerEndpoint == 'http://foo.com'
            assert (it[0] as SubmitContainerTokenRequest).workflowId == '1234'
        }

        cleanup:
        SysEnv.pop()
    }

    // --== launch web server ==--
    @Shared
    def CONFIG_RESP = [
            '/foo.json': new ContainerConfig(entrypoint: ['entry.sh']),
            '/bar.json': new ContainerConfig(layers: [ new ContainerLayer(location: 'http://somewhere.com', gzipDigest:'sha256:1234', gzipSize:100, tarDigest:'sha256:5678') ]),
            'combined': new ContainerConfig(entrypoint: ['entry.sh'], layers: [ new ContainerLayer(location: 'http://somewhere.com', gzipDigest:'sha256:1234', gzipSize:100, tarDigest:'sha256:5678') ]),
    ]

    def 'should resolve container config' () {
        given:
        HttpHandler handler = { HttpExchange exchange ->
            final key = exchange.requestURI.path
            log.debug "Request key: '$key' [$key.class.name] contains: ${CONFIG_RESP.containsKey(key)}"
            final resp = CONFIG_RESP.get(key)
            if( resp ) {
                final body = JsonOutput.toJson(resp)
                exchange.getResponseHeaders().add("Content-Type", "text/json")
                exchange.sendResponseHeaders(200, body.size())
                exchange.getResponseBody() << body
                exchange.getResponseBody().close()
            }
            else {
                log.error("Cannot find response for path '$exchange.requestURI.path'")
                exchange.sendResponseHeaders(404, 0)
            }

        }

        HttpServer server = HttpServer.create(new InetSocketAddress(9901), 0);
        server.createContext("/", handler);
        server.start()

        def session = Mock(Session) {getConfig() >> CONFIG }
        def client = new WaveClient(session)

        expect:
        client.resolveContainerConfig() == EXPECTED

        cleanup:
        server?.stop(0)

        where:
        CONFIG                                      | EXPECTED
        [:]                                         | null
        [ wave:[containerConfigUrl: 'http://localhost:9901/foo.json']] \
                                                    | CONFIG_RESP.get('/foo.json')
        and:
        [ wave:[containerConfigUrl: 'http://localhost:9901/foo.json'], \
          fusion: [enabled: true, containerConfigUrl: 'http://localhost:9901/bar.json']]\
                                                    | CONFIG_RESP.get('combined')
    }

    @Unroll
    def 'should get fusion default url' () {
        given:
        def sess = Mock(Session) {getConfig() >> [:] }
        and:
        def wave = Spy(new WaveClient(sess))

        expect:
        wave.defaultFusionUrl(ARCH).toURI().toString() == EXPECTED
        
        where:
        ARCH                | EXPECTED
        'linux/amd64'       | 'https://fusionfs.seqera.io/releases/v2.2-amd64.json'
        'linux/x86_64'      | 'https://fusionfs.seqera.io/releases/v2.2-amd64.json'
        'arm64'             | 'https://fusionfs.seqera.io/releases/v2.2-arm64.json'
        'linux/arm64'       | 'https://fusionfs.seqera.io/releases/v2.2-arm64.json'
        'linux/arm64/v8'    | 'https://fusionfs.seqera.io/releases/v2.2-arm64.json'
    }

    def 'should check is local conda file' () {
        expect:
        WaveClient.isCondaLocalFile(CONTENT) == EXPECTED

        where:
        CONTENT             | EXPECTED
        'foo'               | false
        'foo.yml'           | true
        'foo.txt'           | true
        'foo\nbar.yml'      | false
        'http://foo.com'    | false
    }

    def 'should check is remote conda file' () {
        expect:
        WaveClient.isCondaRemoteFile(CONTENT) == EXPECTED

        where:
        CONTENT             | EXPECTED
        'foo'               | false
        'foo.yml'           | false
        'foo.txt'           | false
        'foo\nbar.yml'      | false
        'http://foo.com'    | true
        'https://foo.com'   | true
    }

    def 'should retry http request' () {

        given:
        int requestCount=0
        HttpHandler handler = { HttpExchange exchange ->
            if( ++requestCount<3 ) {
                exchange.getResponseHeaders().add("Content-Type", "text/plain")
                exchange.sendResponseHeaders(503, 0)
                exchange.getResponseBody().close()
            }
            else {
                def body = 'Hello world!'
                exchange.getResponseHeaders().add("Content-Type", "text/plain")
                exchange.sendResponseHeaders(200, body.size())
                exchange.getResponseBody().write(body.bytes)
                exchange.getResponseBody().close()
            }
        }

        HttpServer server = HttpServer.create(new InetSocketAddress(9901), 0);
        server.createContext("/", handler);
        server.start()

        def session = Mock(Session) {getConfig() >> [:] }
        def client = new WaveClient(session)

        when:
        def request = HttpRequest.newBuilder().uri(new URI('http://localhost:9901/foo.txt')).build()
        def response = client.httpSend(request)
        then:
        response.statusCode() == 200
        response.body() == 'Hello world!'
        and:
        requestCount == 3
        
        cleanup:
        server?.stop(0)

    }


}
