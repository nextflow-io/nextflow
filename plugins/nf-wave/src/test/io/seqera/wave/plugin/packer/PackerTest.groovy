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

package io.seqera.wave.plugin.packer

import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.attribute.FileTime

import nextflow.extension.FilesEx
import nextflow.file.FileHelper
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class PackerTest extends Specification implements TarHelper {

    def 'should tar bundle' () {
        given:
        def folder = Files.createTempDirectory('test')
        and:
        def result = folder.resolve('result')
        def result2 = folder.resolve('result2')
        and:
        def rootPath = folder.resolve('bundle'); rootPath.mkdir()
        rootPath.resolve('main.nf').text = "I'm the main file"
        rootPath.resolve('this/that').mkdirs()
        Files.write(rootPath.resolve('this/hola.txt'), "Hola".bytes)
        Files.write(rootPath.resolve('this/hello.txt'), "Hello".bytes)
        Files.write(rootPath.resolve('this/that/ciao.txt'), "Ciao".bytes)
        and:
        def files = new ArrayList<Path>()
        FileHelper.visitFiles([type:'any'], rootPath, '**', {
            final mode = it.isDirectory() ? 0700 : 0600
            FilesEx.setPermissionsMode(it, mode)
            //
            files.add(it)
        })
        and:
        def packer = new Packer()

        when:
        def buffer = new ByteArrayOutputStream()
        packer.makeTar(rootPath, files, buffer)
        and:
        untar( new ByteArrayInputStream(buffer.toByteArray()), result )
        then:
        result.resolve('main.nf').text == rootPath.resolve('main.nf').text
        result.resolve('this/hola.txt').text == rootPath.resolve('this/hola.txt').text
        result.resolve('this/hello.txt').text == rootPath.resolve('this/hello.txt').text
        result.resolve('this/that/ciao.txt').text == rootPath.resolve('this/that/ciao.txt').text
        and:
        result.resolve('main.nf').getPermissionsMode() == 0600
        result.resolve('this/hola.txt').getPermissionsMode() == 0600
        result.resolve('this/that').getPermissionsMode() == 0700
        and:
        Files.getLastModifiedTime(result.resolve('main.nf')) == FileTime.fromMillis(0)

        when:
        def layer = packer.layer(rootPath, files)
        then:
        layer.tarDigest == 'sha256:f556b94e9b6f5f72b86e44833614b465df9f65cb4210e3f4416292dca1618360'
        layer.gzipDigest == 'sha256:e58685a82452a11faa926843e7861c94bdb93e2c8f098b5c5354ec9b6fee2b68'
        layer.gzipSize == 251
        and:
        def gzip = layer.location.replace('data:','').decodeBase64()
        def tar = uncompress(gzip)
        untar( new ByteArrayInputStream(tar), result2)
        and:
        result2.resolve('main.nf').text == rootPath.resolve('main.nf').text
        result2.resolve('this/hola.txt').text == rootPath.resolve('this/hola.txt').text
        result2.resolve('this/hello.txt').text == rootPath.resolve('this/hello.txt').text
        result2.resolve('this/that/ciao.txt').text == rootPath.resolve('this/that/ciao.txt').text
        and:
        Files.getLastModifiedTime(result2.resolve('main.nf')) == FileTime.fromMillis(0)

        cleanup:
        folder?.deleteDir()
    }


    def 'should tar bundle and preserve timestamps' () {
        given:
        def LAST_MODIFIED = FileTime.fromMillis(1_000_000_000_000)
        def folder = Files.createTempDirectory('test')
        and:
        def result = folder.resolve('result')
        def result2 = folder.resolve('result2')
        and:
        def rootPath = folder.resolve('bundle'); rootPath.mkdir()
        rootPath.resolve('main.nf').text = "I'm the main file"
        rootPath.resolve('this/that').mkdirs()
        Files.write(rootPath.resolve('this/hola.txt'), "Hola".bytes)
        Files.write(rootPath.resolve('this/hello.txt'), "Hello".bytes)
        Files.write(rootPath.resolve('this/that/ciao.txt'), "Ciao".bytes)
        and:
        def files = new ArrayList<Path>()
        FileHelper.visitFiles([type:'any'], rootPath, '**', {
            Files.setLastModifiedTime(it, LAST_MODIFIED)
            final mode = it.isDirectory() ? 0700 : 0600
            FilesEx.setPermissionsMode(it, mode)
            //
            files.add(it)
        })
        and:
        def packer = new Packer(preserveFileTimestamp: true)

        when:
        def buffer = new ByteArrayOutputStream()
        packer.makeTar(rootPath, files, buffer)
        and:
        untar( new ByteArrayInputStream(buffer.toByteArray()), result )
        then:
        result.resolve('main.nf').text == rootPath.resolve('main.nf').text
        result.resolve('this/hola.txt').text == rootPath.resolve('this/hola.txt').text
        result.resolve('this/hello.txt').text == rootPath.resolve('this/hello.txt').text
        result.resolve('this/that/ciao.txt').text == rootPath.resolve('this/that/ciao.txt').text
        and:
        result.resolve('main.nf').getPermissionsMode() == 0600
        result.resolve('this/hola.txt').getPermissionsMode() == 0600
        result.resolve('this/that').getPermissionsMode() == 0700
        and:
        Files.getLastModifiedTime(result.resolve('main.nf')) == LAST_MODIFIED

        when:
        def layer = packer.layer(rootPath, files)
        then:
        layer.tarDigest == 'sha256:81200f6ad32793567d8070375dc51312a1711fedf6a1c6f5e4a97fa3014f3491'
        layer.gzipDigest == 'sha256:09a2deca4293245909223db505cf69affa1a8ff8acb745fe3cad38bc0b719110'
        layer.gzipSize == 254
        and:
        def gzip = layer.location.replace('data:','').decodeBase64()
        def tar = uncompress(gzip)
        untar( new ByteArrayInputStream(tar), result2)
        and:
        result2.resolve('main.nf').text == rootPath.resolve('main.nf').text
        result2.resolve('this/hola.txt').text == rootPath.resolve('this/hola.txt').text
        result2.resolve('this/hello.txt').text == rootPath.resolve('this/hello.txt').text
        result2.resolve('this/that/ciao.txt').text == rootPath.resolve('this/that/ciao.txt').text

        cleanup:
        folder?.deleteDir()
    }

}
