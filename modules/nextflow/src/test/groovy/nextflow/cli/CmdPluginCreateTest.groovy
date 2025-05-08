/*
 * Copyright 2013-2025, Seqera Labs
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

import java.nio.file.Files
import java.nio.file.Path

import nextflow.plugin.Plugins
import spock.lang.IgnoreIf
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CmdPluginCreateTest extends Specification {

    def cleanup() {
        Plugins.stop()
    }

    @IgnoreIf({System.getenv('NXF_SMOKE')})
    def 'should clone and create a plugin project' () {
        given:
        def folder = Files.createTempDirectory('test')
        and:
        def args = [
            'create',
            'hello world plugin',
            'foo',
            folder.toAbsolutePath().toString() + '/hello']

        when:
        def cmd = new CmdPlugin(args: args)
        and:
        cmd.run()

        then:
        Files.exists(folder.resolve('hello'))
        Files.exists(folder.resolve('hello/src/main/groovy/foo/plugin/HelloWorldPlugin.groovy'))
        Files.exists(folder.resolve('hello/src/main/groovy/foo/plugin/HelloWorldObserver.groovy'))
        Files.exists(folder.resolve('hello/src/main/groovy/foo/plugin/HelloWorldFactory.groovy'))
        and:
        Files.exists(folder.resolve('hello/src/test/groovy/foo/plugin/HelloWorldObserverTest.groovy'))
        and:
        Path.of(folder.resolve('hello/settings.gradle').toUri()).text.contains("rootProject.name = 'hello-world-plugin'")
        Path.of(folder.resolve('hello/build.gradle').toUri()).text.contains("provider = 'foo'")
        Path.of(folder.resolve('hello/build.gradle').toUri()).text.contains("className = 'foo.plugin.HelloWorldPlugin'")
        and:
        Path.of(folder.resolve('hello/README.md').toUri()).text.contains("# hello-world-plugin plugin")
        Path.of(folder.resolve('hello/README.md').toUri()).text.contains("nextflow run hello -plugins hello-world-plugin@0.1.0")

        cleanup:
        folder?.deleteDir()
    }

}
