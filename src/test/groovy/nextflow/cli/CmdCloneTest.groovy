/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
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
        cmd.args = ['nextflow-io/hello', dir.toFile().toString(), ]

        when:
        cmd.run()

        then:
        dir.resolve('README.md').exists()

        cleanup:
        dir?.deleteDir()

    }

}
