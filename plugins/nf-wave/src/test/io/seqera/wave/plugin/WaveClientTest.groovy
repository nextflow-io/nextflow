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

import groovy.transform.CompileStatic
import io.seqera.wave.plugin.WaveClient
import nextflow.Session
import nextflow.extension.FilesEx
import nextflow.file.FileHelper
import nextflow.script.bundle.ModuleBundle
import org.apache.commons.compress.archivers.ArchiveStreamFactory
import org.apache.commons.compress.archivers.tar.TarArchiveEntry
import org.apache.commons.compress.archivers.tar.TarArchiveInputStream
import org.apache.commons.compress.compressors.gzip.GzipCompressorInputStream
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
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
        def result = folder.resolve('result')
        def result2 = folder.resolve('result2')
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
        def bundle = ModuleBundle.scan(bundlePath)

        when:
        def wave = new WaveClient(sess)
        def buffer = new ByteArrayOutputStream()
        wave.makeTar(bundle, buffer)
        and:
        untar( new ByteArrayInputStream(buffer.toByteArray()), result )
        then:
        result.resolve('main.nf').text == bundlePath.resolve('main.nf').text
        result.resolve('this/hola.txt').text == bundlePath.resolve('this/hola.txt').text
        result.resolve('this/hello.txt').text == bundlePath.resolve('this/hello.txt').text
        result.resolve('this/that/ciao.txt').text == bundlePath.resolve('this/that/ciao.txt').text
        and:
        result.resolve('main.nf').getPermissionsMode() == 0600
        result.resolve('this/hola.txt').getPermissionsMode() == 0600
        result.resolve('this/that').getPermissionsMode() == 0700
        and:
        Files.getLastModifiedTime(result.resolve('main.nf')) == LAST_MODIFIED
        
        when:
        def layer = wave.layer(bundle)
        then:
        layer.tarDigest == 'sha256:81200f6ad32793567d8070375dc51312a1711fedf6a1c6f5e4a97fa3014f3491'
        layer.gzipDigest == 'sha256:09a2deca4293245909223db505cf69affa1a8ff8acb745fe3cad38bc0b719110'
        layer.gzipSize == 254
        and:
        def gzip = layer.location.replace('data:','').decodeBase64()
        def tar = uncompress(gzip)
        untar( new ByteArrayInputStream(tar), result2)
        and:
        result2.resolve('main.nf').text == bundlePath.resolve('main.nf').text
        result2.resolve('this/hola.txt').text == bundlePath.resolve('this/hola.txt').text
        result2.resolve('this/hello.txt').text == bundlePath.resolve('this/hello.txt').text
        result2.resolve('this/that/ciao.txt').text == bundlePath.resolve('this/that/ciao.txt').text

        cleanup:
        folder?.deleteDir()
    }

}
