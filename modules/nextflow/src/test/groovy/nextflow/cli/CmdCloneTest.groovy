/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

import spock.lang.Requires
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CmdCloneTest extends Specification {

    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def testClone() {

        given:
        def accessToken = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def dir = Files.createTempDirectory('test')
        def cmd = new CmdClone(hubUser: accessToken)
        cmd.args = ['nextflow-io/hello', dir.toFile().toString()]

        when:
        cmd.run()

        then:
        dir.resolve('README.md').exists()

        cleanup:
        dir?.deleteDir()

    }

}
