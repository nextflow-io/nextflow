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

package io.seqera.wave.plugin.cli

import java.nio.file.Files
import java.nio.file.Path

import groovy.json.JsonSlurper
import io.seqera.wave.plugin.packer.TarHelper
import nextflow.extension.FilesEx
import nextflow.file.FileHelper
import org.junit.Rule
import spock.lang.Shared
import spock.lang.Specification
import spock.lang.TempDir
import spock.lang.Unroll
import test.OutputCapture
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class WaveCmdEntryTest extends Specification implements TarHelper {

    @Shared
    @TempDir
    Path folder

    @Rule
    OutputCapture capture = new OutputCapture()

    def 'should create create pack' () {
        given:
        def rootPath = folder.resolve('bundle'); rootPath.mkdir()
        def untarPath = folder.resolve('untar'); untarPath.mkdir()
        rootPath.resolve('main.nf').text = "I'm the main file"
        rootPath.resolve('this/that').mkdirs()
        Files.write(rootPath.resolve('this/hola.txt'), "Hola".bytes)
        Files.write(rootPath.resolve('this/hello.txt'), "Hello".bytes)
        Files.write(rootPath.resolve('this/that/ciao.txt'), "Ciao".bytes)
        and:
        FileHelper.visitFiles([type:'any'], rootPath, '**', {
            final mode = it.isDirectory() ? 0700 : 0600
            FilesEx.setPermissionsMode(it, mode)
        })
        and:
        def cmd = new WaveCmdEntry()

        when:
        def result = cmd.packContainer( [rootPath.toString()] )
        then:
        def gzipFile = folder.resolve('bundle.tar.gz')
        gzipFile.exists()
        and:
        def json = new JsonSlurper().parseText(result)
        and:
        json.layers[0].gzipSize == Files.size(gzipFile)
        json.layers[0].location == gzipFile.toUri().toString()
        json.layers[0].tarDigest == 'sha256:f556b94e9b6f5f72b86e44833614b465df9f65cb4210e3f4416292dca1618360'
        json.layers[0].gzipDigest == 'sha256:e58685a82452a11faa926843e7861c94bdb93e2c8f098b5c5354ec9b6fee2b68'
        
        when:
        def tar = uncompress(Files.readAllBytes(gzipFile))
        untar( new ByteArrayInputStream(tar), untarPath )
        then:
        untarPath.resolve('main.nf').text == rootPath.resolve('main.nf').text
        untarPath.resolve('this/hola.txt').text == rootPath.resolve('this/hola.txt').text
        untarPath.resolve('this/hello.txt').text == rootPath.resolve('this/hello.txt').text
        untarPath.resolve('this/that/ciao.txt').text == rootPath.resolve('this/that/ciao.txt').text
    }

    @Unroll
    def 'should find base name' () {
        expect:
        WaveCmdEntry.baseName(PATH) == EXPECTED

        where:
        PATH                        | EXPECTED
        null                        | null
        '/some/name'                | 'name'
        'http://some/name.tar'      | 'name'
        '/some/name.tar.gz'         | 'name'
        '/some/name.tar.gzip'       | 'name'
        'http://some/name.tar.gz'   | 'name'
        'http://some/'              | null
        'http://some'               | null
    }
}
