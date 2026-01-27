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

import nextflow.exception.AbortOperationException
import nextflow.plugin.Plugins
import spock.lang.Specification

/**
 * Tests for CmdPush command
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */

class CmdPushTest extends Specification {

    def cleanup() {
        Plugins.stop()
    }

    def 'should fail when no repository specified and no git repo exists'() {
        given:
        def tempDir = Files.createTempDirectory('test').toFile()
        def cmd = new CmdPush()
        cmd.rootFolder = tempDir

        when:
        cmd.run()

        then:
        def e = thrown(AbortOperationException)
        e.message.contains('No git repository found')

        cleanup:
        tempDir?.deleteDir()
    }

    def 'should get command name'() {
        given:
        def cmd = new CmdPush()

        expect:
        cmd.getName() == 'push'
    }

    def 'should have default parameter values'() {
        given:
        def cmd = new CmdPush()

        expect:
        cmd.maxSizeMB == 10
        cmd.message == 'Push from nextflow'
    }
}
