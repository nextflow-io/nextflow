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
 */

package nextflow.cli

import static nextflow.scm.AssetManager.REVISION_SUBDIR

import nextflow.plugin.Plugins
import spock.lang.IgnoreIf

import java.nio.file.Files

import spock.lang.Requires
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@IgnoreIf({System.getenv('NXF_SMOKE')})
class CmdPullTest extends Specification {

    def cleanup() {
        Plugins.stop()
    }
    
    @Requires({ System.getenv('NXF_GITHUB_ACCESS_TOKEN') })
    def 'should pull the github repository in the local folder'() {

        given:
        def accessToken = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def dir = Files.createTempDirectory('test')
        def cmd = new CmdPull(args: ['nextflow-io/hello'], root: dir.toFile(), revision: '7588c46ffefb4e3c06d4ab32c745c4d5e56cdad8', hubUser: accessToken)

        when:
        cmd.run()
        then:
        dir.resolve('nextflow-io/hello/' + REVISION_SUBDIR + '/' + '7588c46ffefb4e3c06d4ab32c745c4d5e56cdad8' + '/.git').exists()
        dir.resolve('nextflow-io/hello/' + REVISION_SUBDIR + '/' + '7588c46ffefb4e3c06d4ab32c745c4d5e56cdad8' + '/README.md').exists()

        cleanup:
        dir?.deleteDir()
        Plugins.stop()
    }

}
