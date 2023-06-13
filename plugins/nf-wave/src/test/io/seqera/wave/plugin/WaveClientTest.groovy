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
        and:
        req.fingerprint == 'bd2cb4b32df41f2d290ce2366609f2ad'
        req.timestamp instanceof String
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
        def assets = client.resolveAssets(task, IMAGE)
        then:
        assets.containerImage == IMAGE
        !assets.moduleResources
        !assets.dockerFileContent
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
        def assets = client.resolveAssets(task, IMAGE)
        then:
        assets.containerImage == IMAGE
        assets.containerPlatform == 'linux/arm64'
        !assets.moduleResources
        !assets.dockerFileContent
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
        def assets = client.resolveAssets(task, IMAGE)
        then:
        assets.containerImage == IMAGE
        assets.moduleResources == BUNDLE
        !assets.dockerFileContent
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
        def assets = client.resolveAssets(task, IMAGE)
        then:
        client.resolveContainerConfig(ARCH) >> CONTAINER_CONFIG
        and:
        assets.containerImage == IMAGE
        assets.moduleResources == BUNDLE
        assets.containerConfig == CONTAINER_CONFIG
        and:
        !assets.dockerFileContent
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
        def assets = client.resolveAssets(task, null)
        then:
        assets.dockerFileContent == 'FROM foo\nRUN this/that'
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
        def task = Mock(TaskRun) {getConfig() >> [conda:'salmon=1.2.3'] }
        and:
        def client = new WaveClient(session)

        when:
        def assets = client.resolveAssets(task, null)
        then:
        assets.dockerFileContent == '''\
                FROM mambaorg/micromamba:1.4.2
                RUN \\
                    micromamba install -y -n base -c conda-forge -c defaults \\
                    salmon=1.2.3 \\
                    && micromamba clean -a -y
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
        def task = Mock(TaskRun) {getConfig() >> [spack:'salmon@1.2.3', arch:'amd64'] }
        and:
        def client = new WaveClient(session)

        when:
        def assets = client.resolveAssets(task, null)
        then:
        assets.dockerFileContent == '''\
# Builder image
FROM spack/ubuntu-jammy:v0.19.2 as builder

RUN mkdir -p /opt/spack-env \\
&&  spack env create -d /opt/spack-env \\
&&  sed -e 's;compilers:;compilers::;' \\
         -e 's;^ *flags: *{};    flags:\\n      cflags: -O3\\n      cxxflags: -O3\\n      fflags: -O3;' \\
         /root/.spack/linux/compilers.yaml > /opt/spack-env/compilers.yaml \\
&&  sed -i '/^spack:/a\\  include: [/opt/spack-env/compilers.yaml]' /opt/spack-env/spack.yaml \\
&& cd /opt/spack-env && spack env activate . \\
&& spack add salmon@1.2.3 \\
&& spack config add config:install_tree:/opt/software \\
&& spack config add concretizer:unify:true \\
&& spack config add concretizer:reuse:false \\
&& spack config add packages:all:target:[x86_64] \\
&& echo -e "\\
  view: /opt/view \\n\\
" >> /opt/spack-env/spack.yaml

# Install packages, clean afterwards, finally strip binaries
RUN cd /opt/spack-env && spack env activate . \\
&& spack concretize -f \\
&& spack install --fail-fast && spack gc -y \\
&& find -L /opt/._view/* -type f -exec readlink -f '{}' \\; | \\
    xargs file -i | \\
    grep 'charset=binary' | \\
    grep 'x-executable\\|x-archive\\|x-sharedlib' | \\
    awk -F: '{print \$1}' | xargs strip -s

RUN cd /opt/spack-env && \\
    spack env activate --sh -d . >> /opt/spack-env/z10_spack_environment.sh && \\
    original_view=\$( cd /opt ; ls -1d ._view/* ) && \\
    sed -i "s;/view/;/\$original_view/;" /opt/spack-env/z10_spack_environment.sh && \\
    echo "# Needed for Perl applications" >>/opt/spack-env/z10_spack_environment.sh && \\
    echo "export PERL5LIB=\$(eval ls -d /opt/._view/*/lib/5.*):\$PERL5LIB" >>/opt/spack-env/z10_spack_environment.sh && \\
    rm -rf /opt/view

# Runner image
FROM ubuntu:22.04

COPY --from=builder /opt/spack-env /opt/spack-env
COPY --from=builder /opt/software /opt/software
COPY --from=builder /opt/._view /opt/._view

# Near OS-agnostic package addition
RUN ( apt update -y && apt install -y procps libgomp1  && rm -rf /var/lib/apt/lists/* ) || \\
    ( yum install -y procps libgomp  && yum clean all && rm -rf /var/cache/yum ) || \\
    ( zypper ref && zypper install -y procps libgomp1  && zypper clean -a ) || \\
    ( apk update && apk add --no-cache procps libgomp bash  && rm -rf /var/cache/apk )

# Entrypoint for Singularity
RUN mkdir -p /.singularity.d/env && \\
    cp -p /opt/spack-env/z10_spack_environment.sh /.singularity.d/env/91-environment.sh
# Entrypoint for Docker
RUN echo "#!/usr/bin/env bash\\n\\nset -ef -o pipefail\\nsource /opt/spack-env/z10_spack_environment.sh\\nexec \\"\\\$@\\"" \\
    >/opt/spack-env/spack_docker_entrypoint.sh && chmod a+x /opt/spack-env/spack_docker_entrypoint.sh


ENTRYPOINT [ "/opt/spack-env/spack_docker_entrypoint.sh" ]
CMD [ "/bin/bash" ]
'''//.stripIndent()
        and:
        !assets.moduleResources
        !assets.containerImage
        !assets.containerConfig
        !assets.condaFile
        !assets.spackFile
        !assets.projectResources
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
        def assets = client.resolveAssets(task, null)
        then:
        assets.dockerFileContent == '''\
                FROM mambaorg/micromamba:1.4.2
                COPY --chown=$MAMBA_USER:$MAMBA_USER conda.yml /tmp/conda.yml
                RUN micromamba install -y -n base -f /tmp/conda.yml \\
                    && micromamba clean -a -y
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
        def assets = client.resolveAssets(task, null)
        then:
        assets.dockerFileContent == '''\
# Builder image
FROM spack/ubuntu-jammy:v0.19.2 as builder
COPY spack.yaml /tmp/spack.yaml

RUN mkdir -p /opt/spack-env \\
&&  sed -e 's;compilers:;compilers::;' \\
         -e 's;^ *flags: *{};    flags:\\n      cflags: -O3\\n      cxxflags: -O3\\n      fflags: -O3;' \\
         /root/.spack/linux/compilers.yaml > /opt/spack-env/compilers.yaml \\
&&  sed '/^spack:/a\\  include: [/opt/spack-env/compilers.yaml]' /tmp/spack.yaml > /opt/spack-env/spack.yaml \\
&& cd /opt/spack-env && spack env activate . \\
&& spack config add config:install_tree:/opt/software \\
&& spack config add concretizer:unify:true \\
&& spack config add concretizer:reuse:false \\
&& spack config add packages:all:target:[x86_64] \\
&& echo -e "\\
  view: /opt/view \\n\\
" >> /opt/spack-env/spack.yaml

# Install packages, clean afterwards, finally strip binaries
RUN cd /opt/spack-env && spack env activate . \\
&& spack concretize -f \\
&& spack install --fail-fast && spack gc -y \\
&& find -L /opt/._view/* -type f -exec readlink -f '{}' \\; | \\
    xargs file -i | \\
    grep 'charset=binary' | \\
    grep 'x-executable\\|x-archive\\|x-sharedlib' | \\
    awk -F: '{print \$1}' | xargs strip -s

RUN cd /opt/spack-env && \\
    spack env activate --sh -d . >> /opt/spack-env/z10_spack_environment.sh && \\
    original_view=\$( cd /opt ; ls -1d ._view/* ) && \\
    sed -i "s;/view/;/\$original_view/;" /opt/spack-env/z10_spack_environment.sh && \\
    echo "# Needed for Perl applications" >>/opt/spack-env/z10_spack_environment.sh && \\
    echo "export PERL5LIB=\$(eval ls -d /opt/._view/*/lib/5.*):\$PERL5LIB" >>/opt/spack-env/z10_spack_environment.sh && \\
    rm -rf /opt/view

# Runner image
FROM ubuntu:22.04

COPY --from=builder /opt/spack-env /opt/spack-env
COPY --from=builder /opt/software /opt/software
COPY --from=builder /opt/._view /opt/._view

# Near OS-agnostic package addition
RUN ( apt update -y && apt install -y procps libgomp1  && rm -rf /var/lib/apt/lists/* ) || \\
    ( yum install -y procps libgomp  && yum clean all && rm -rf /var/cache/yum ) || \\
    ( zypper ref && zypper install -y procps libgomp1  && zypper clean -a ) || \\
    ( apk update && apk add --no-cache procps libgomp bash  && rm -rf /var/cache/apk )

# Entrypoint for Singularity
RUN mkdir -p /.singularity.d/env && \\
    cp -p /opt/spack-env/z10_spack_environment.sh /.singularity.d/env/91-environment.sh
# Entrypoint for Docker
RUN echo "#!/usr/bin/env bash\\n\\nset -ef -o pipefail\\nsource /opt/spack-env/z10_spack_environment.sh\\nexec \\"\\\$@\\"" \\
    >/opt/spack-env/spack_docker_entrypoint.sh && chmod a+x /opt/spack-env/spack_docker_entrypoint.sh


ENTRYPOINT [ "/opt/spack-env/spack_docker_entrypoint.sh" ]
CMD [ "/bin/bash" ]
'''//.stripIndent()
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
        def assets = wave.resolveAssets(task, 'image:latest')
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

    @Unroll
    def 'should find wave token' () {
        given:
        def sess = Mock(Session) {getConfig() >> [wave:[endpoint: 'http://foo.com']] }
        and:
        def wave = Spy(new WaveClient(sess))

        expect:
        wave.getWaveToken(CONTAINER) == EXPECTED

        where:
        EXPECTED            | CONTAINER
        null                | null
        null                | 'ubunutu:latest'
        null                | 'xyz.com/wt/3aec54700cff/wave-build-repo:rnaseq-nf_v1.0'
        and:
        '3aec54700cff'      | 'foo.com/wt/3aec54700cff/wave-build-repo:rnaseq-nf_v1.0'
    }

    def 'should convert json to describe container response'() {
        given:
        def JSON = '''
            {
              "token": "3aec54700cff",
              "expiration": "2023-03-02T18:07:50.488226285Z",
              "request": {
                "user": {
                  "id": 8083,
                  "userName": "pditommaso",
                  "email": "pditommaso@me.com"
                },
                "workspaceId": 88265370860066,
                "containerImage": "1234567890.dkr.ecr.us-west-2.amazonaws.com/wave-build-repo:rnaseq-nf_v1.0",
                "containerConfig": {
                  "entrypoint": [
                    "/opt/fusion/entrypoint.sh"
                  ],
                  "layers": [
                    {
                      "location": "data:DATA+OMITTED",
                      "gzipDigest": "sha256:dc8dd4ebf839869abb81d35fe3f265de9f3ac7b9b285e274c6b92072b02a84ec",
                      "gzipSize": 202,
                      "tarDigest": "sha256:dc4d652cd223da5bca40d08890686c4198769fb7bfc09de2ed3c3c77dead4bf9"
                    },
                    {
                      "location": "https://fusionfs.seqera.io/releases/pkg/0/6/4/fusionfs-amd64.tar.gz",
                      "gzipDigest": "sha256:c55640ae3284715e5c7a1c1f6c6ec2de77a881a08f5a9c46f077ecd0379e8477",
                      "gzipSize": 6191418,
                      "tarDigest": "sha256:e24642d65d5b21987666cf1ce4ba007ecadedbcefae9601669ab43a566682aa6"
                    }
                  ]
                },
                "towerEndpoint": "https://api.tower.nf",
                "fingerprint": "779855a0ffc582ef3170f7dab8829465",
                "timestamp": "2023-03-01T20:07:49.933811174Z",
                "zoneId": "Z",
                "ipAddress": "54.190.237.226"
              },
              "build": {
                "buildRepository": "1234567890.dkr.ecr.us-west-2.amazonaws.com/wave-build-repo",
                "cacheRepository": "1234567890.dkr.ecr.us-west-2.amazonaws.com/wave-cache-repo"
              },
              "source": {
                "image": "1234567890.dkr.ecr.us-west-2.amazonaws.com/wave-build-repo:rnaseq-nf_v1.0",
                "digest": "sha256:d6f56ed0eae171fabd324bf582dd5c49c6462662c80a7e69632c57043b6af143"
              },
              "wave": {
                "image": "wave.seqera.io/wt/3aec54700cff/wave-build-repo:rnaseq-nf_v1.0",
                "digest": "sha256:d8f4f9aa77b4d1941b50a050ed71473a0e04720f38a12d497557c39a25398830"
              }
            }
            '''
        and:
        def sess = Mock(Session) {getConfig() >> [:] }
        def wave = Spy(new WaveClient(sess))

        when:
        def resp = wave.jsonToDescribeContainerResponse(JSON)
        then:
        resp.token == '3aec54700cff'
        and:
        resp.wave.image == 'wave.seqera.io/wt/3aec54700cff/wave-build-repo:rnaseq-nf_v1.0'
        resp.wave.digest == 'sha256:d8f4f9aa77b4d1941b50a050ed71473a0e04720f38a12d497557c39a25398830'
        and:
        resp.request.user.id == 8083
        resp.request.user.userName == 'pditommaso'

    }

    def 'should resolve wave container' () {
        given:
        def RESP1 = '''
            {
              "token": "3aec54700cff",
              "source": {
                "image": "docker.io/library/ubuntu:latest",
                "digest": "sha256:d6f56ed0eae171fabd324bf582dd5c49c6462662c80a7e69632c57043b6af143"
              },
              "wave": {
                "image": "wave.seqera.io/wt/3aec54700cff/library/ubuntu:latest",
                "digest": "sha256:d8f4f9aa77b4d1941b50a050ed71473a0e04720f38a12d497557c39a25398830"
              }
            }
        '''
        and:
        def sess = Mock(Session) {getConfig() >> [:] }
        def wave = Spy(new WaveClient(sess))
        
        when:
        def result = wave.resolveSourceContainer('ubuntu')
        then:
        0 * wave.fetchContainerInfo(_) >> null
        result == 'ubuntu'

        when:
        result = wave.resolveSourceContainer('wave.seqera.io/wt/3aec54700cff/library/ubuntu:latest')
        then:
        1 * wave.fetchContainerInfo('3aec54700cff') >> RESP1
        result == 'wave.seqera.io/wt/3aec54700cff/library/ubuntu@sha256:d8f4f9aa77b4d1941b50a050ed71473a0e04720f38a12d497557c39a25398830'
    }

    def 'should return source container' () {
        given:
        def RESP = '''
            {
              "token": "3aec54700cff",
              "source": {
                "image": "docker.io/library/ubuntu:latest",
                "digest": "sha256:d6f56ed0eae171fabd324bf582dd5c49c6462662c80a7e69632c57043b6af143"
              },
              "wave": {
                "image": "wave.seqera.io/wt/3aec54700cff/library/ubuntu:latest",
                "digest": "sha256:d6f56ed0eae171fabd324bf582dd5c49c6462662c80a7e69632c57043b6af143"
              }
            }
        '''
        and:
        def sess = Mock(Session) {getConfig() >> [:] }
        def wave = Spy(new WaveClient(sess))

        when:
        def result = wave.resolveSourceContainer('wave.seqera.io/wt/3aec54700cff/library/ubuntu:latest')
        then:
        1 * wave.fetchContainerInfo('3aec54700cff') >> RESP
        result == 'docker.io/library/ubuntu@sha256:d6f56ed0eae171fabd324bf582dd5c49c6462662c80a7e69632c57043b6af143'
    }
}
