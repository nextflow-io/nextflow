/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

package nextflow.cli

import spock.lang.Requires
import spock.lang.Shared
import spock.lang.Specification

import java.nio.file.Files
import java.nio.file.Path

import groovy.json.JsonSlurper
import nextflow.scm.AssetManager
import org.yaml.snakeyaml.Yaml
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
class CmdInfoTest extends Specification {

    @Shared Path tempDir

    def setupSpec() {
        tempDir = Files.createTempDirectory('test')
        AssetManager.root = tempDir.toFile()
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def manager = new AssetManager().build('nextflow-io/hello', [providers: [github: [auth: token]]])
        // download the project
        manager.download()
    }

    def cleanupSpec() {
        tempDir?.deleteDir()
    }

    def 'should print project info' () {

        given:
        def buffer = new ByteArrayOutputStream()

        when:
        new CmdInfo(args: ['hello'], out: new PrintStream(buffer)).run()
        def screen = buffer.toString()

        then:
        screen.contains(" project name: nextflow-io/hello")
        screen.contains(" repository  : https://github.com/nextflow-io/hello")
        screen.contains(" local path  : $tempDir/nextflow-io/hello" )
        screen.contains(" main script : main.nf")
        screen.contains(" revisions   : ")
        screen.contains(" * master (default)")
    }

    def 'should print json info' () {

        given:
        def buffer = new ByteArrayOutputStream()
        def cmd = new CmdInfo(args: ['hello'], format: 'json', out: new PrintStream(buffer))

        when:
        cmd.run()
        def screen = buffer.toString()
        def json = (Map)new JsonSlurper().parseText(screen)

        then:
        json.projectName == "nextflow-io/hello"
        json.repository == "https://github.com/nextflow-io/hello"
        json.localPath == "$tempDir/nextflow-io/hello"
        json.manifest.mainScript == 'main.nf'
        json.manifest.defaultBranch == 'master'
        json.revisions.current == 'master'
        json.revisions.master == 'master'
        json.revisions.branches.size()>1
        json.revisions.branches.any { it.name == 'master' }
        json.revisions.tags.size()>1
        json.revisions.tags.any { it.name == 'v1.1' }

    }

    def 'should print yaml info' () {

        given:
        def buffer = new ByteArrayOutputStream()
        def cmd = new CmdInfo(args: ['hello'], format: 'yaml', out: new PrintStream(buffer))

        when:
        cmd.run()
        def screen = buffer.toString()
        def json = (Map)new Yaml().load(screen)

        then:
        json.projectName == "nextflow-io/hello"
        json.repository == "https://github.com/nextflow-io/hello"
        json.localPath == "$tempDir/nextflow-io/hello"
        json.manifest.mainScript == 'main.nf'
        json.manifest.defaultBranch == 'master'
        json.revisions.current == 'master'
        json.revisions.master == 'master'
        json.revisions.branches.size()>1
        json.revisions.branches.any { it.name == 'master' }
        json.revisions.tags.size()>1
        json.revisions.tags.any { it.name == 'v1.1' }

    }
}
