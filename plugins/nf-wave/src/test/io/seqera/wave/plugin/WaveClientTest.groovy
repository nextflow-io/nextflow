/*
 * Copyright 2020-2022, Seqera Labs
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
import nextflow.conda.CondaConfig
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
        def assets = new WaveAssets('my:image', null, MODULE_RES, null, null, null, PROJECT_RES)
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

    def 'should create dockerfile content from conda recipe' () {
        given:
        def session = Mock(Session) { getConfig() >> [:]}
        def RECIPE = 'bwa=0.7.15 salmon=1.1.1'
        when:
        def client = new WaveClient(session)
        then:
        client.condaRecipeToDockerFile(RECIPE) == '''\
                FROM mambaorg/micromamba:1.0.0
                RUN \\
                   micromamba install -y -n base -c conda-forge -c defaults \\
                   bwa=0.7.15 salmon=1.1.1 \\
                   && micromamba clean -a -y
                '''.stripIndent()
    }

    def 'should create dockerfile content with custom config' () {
        given:
        def session = Mock(Session) {
            getConfig() >> [:]
            getCondaConfig() >> new CondaConfig([channels:'foo,bar'], [:])
        }
        def RECIPE = 'bwa=0.7.15 salmon=1.1.1'
        when:
        def client = new WaveClient(session)
        then:
        client.condaRecipeToDockerFile(RECIPE) == '''\
                FROM mambaorg/micromamba:1.0.0
                RUN \\
                   micromamba install -y -n base -c foo -c bar \\
                   bwa=0.7.15 salmon=1.1.1 \\
                   && micromamba clean -a -y
                '''.stripIndent()
    }

    def 'should create dockerfile content with custom channels' () {
        given:
        def CONDA_OPTS = [mambaImage:'my-base:123', commands: ['USER my-user', 'RUN apt-get update -y && apt-get install -y procps']]
        def session = Mock(Session) { getConfig() >> [wave:[build:[conda:CONDA_OPTS]]]}
        def RECIPE = 'bwa=0.7.15 salmon=1.1.1'
        when:
        def client = new WaveClient(session)
        then:
        client.condaRecipeToDockerFile(RECIPE) == '''\
                FROM my-base:123
                RUN \\
                   micromamba install -y -n base -c conda-forge -c defaults \\
                   bwa=0.7.15 salmon=1.1.1 \\
                   && micromamba clean -a -y
                USER my-user
                RUN apt-get update -y && apt-get install -y procps
                '''.stripIndent()
    }

    def 'should create dockerfile content from conda file' () {
        given:
        def session = Mock(Session) { getConfig() >> [:]}
        when:
        def client = new WaveClient(session)
        then:
        client.condaFileToDockerFile()== '''\
                FROM mambaorg/micromamba:1.0.0
                COPY --chown=$MAMBA_USER:$MAMBA_USER conda.yml /tmp/conda.yml
                RUN micromamba install -y -n base -f /tmp/conda.yml && \\
                    micromamba clean -a -y
                '''.stripIndent()

    }

    def 'should create asset with image' () {
        given:
        def session = Mock(Session) { getConfig() >> [:]}
        def task = Mock(TaskRun) { getConfig() >> [:] }
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
        !assets.projectResources
        !assets.containerPlatform
    }

    def 'should create asset with image and platform' () {
        given:
        def session = Mock(Session) { getConfig() >> [wave:[containerPlatform:'linux/amd64']]}
        def task = Mock(TaskRun) { getConfig() >> [:] }
        def IMAGE = 'foo:latest'
        and:
        def client = new WaveClient(session)

        when:
        def assets = client.resolveAssets(task, IMAGE)
        then:
        assets.containerImage == IMAGE
        assets.containerPlatform == 'linux/amd64'
        !assets.moduleResources
        !assets.dockerFileContent
        !assets.containerConfig
        !assets.condaFile
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
        !assets.projectResources
    }

    def 'should create asset with image and bundle and container config' () {
        given:
        def IMAGE = 'foo:latest'
        def BUNDLE = Mock(ResourcesBundle)
        def CONTAINER_CONFIG = new ContainerConfig(entrypoint: ['entry.sh'], layers: [new ContainerLayer(location: 'http://somewhere')])
        and:
        def session = Mock(Session) { getConfig() >> [:]}
        def task = Mock(TaskRun) { getConfig() >> [:]; getModuleBundle() >> BUNDLE }
        and:
        WaveClient client = Spy(WaveClient, constructorArgs:[session])

        when:
        def assets = client.resolveAssets(task, IMAGE)
        then:
        client.resolveContainerConfig() >> CONTAINER_CONFIG
        and:
        assets.containerImage == IMAGE
        assets.moduleResources == BUNDLE
        assets.containerConfig == CONTAINER_CONFIG
        and:
        !assets.dockerFileContent
        !assets.condaFile
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
                    FROM mambaorg/micromamba:1.0.0
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
                    FROM mambaorg/micromamba:1.0.0
                    COPY --chown=$MAMBA_USER:$MAMBA_USER conda.yml /tmp/conda.yml
                    RUN micromamba install -y -n base -f /tmp/conda.yml && \\
                        micromamba clean -a -y
                    '''.stripIndent()
        and:
        assets.condaFile == condaFile
        and:
        !assets.moduleResources
        !assets.containerImage
        !assets.containerConfig
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
        and:
        def task = Mock(TaskRun) {getModuleBundle() >> MODULE_RES; getConfig() >> [:] }
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
        1 * wave.resolveContainerConfig() >> CONTAINER_CONFIG
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
        e.message == "Process 'foo' declares both 'container' and 'conda' directives that conflicts each other"

        when:
        client.checkConflicts([conda:'this', dockerfile:'that'], 'foo')
        then:
        e = thrown(IllegalArgumentException)
        e.message == "Process 'foo' declares both a 'conda' directive and a module bundle dockerfile that conflicts each other"

        when:
        client.checkConflicts([container:'this', dockerfile:'that'], 'foo')
        then:
        e = thrown(IllegalArgumentException)
        e.message == "Process 'foo' declares both a 'container' directive and a module bundle dockerfile that conflicts each other"

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
        def config = CONFIG
        def sess = Mock(Session) {getConfig() >> config }
        and:
        def wave = Spy(new WaveClient(sess))

        expect:
        wave.defaultFusionUrl().toURI().toString() == EXPECTED
        
        where:
        CONFIG                                      | EXPECTED
        [:]                                         | 'https://fusionfs.seqera.io/releases/v0.6.0-amd64.json'
        [wave:[containerPlatform: 'arm64']]         | 'https://fusionfs.seqera.io/releases/v0.6.0-arm64.json'
        [wave:[containerPlatform: 'linux/arm64']]   | 'https://fusionfs.seqera.io/releases/v0.6.0-arm64.json'
        [wave:[containerPlatform: 'linux/arm64/v8']]    | 'https://fusionfs.seqera.io/releases/v0.6.0-arm64.json'
    }
}
