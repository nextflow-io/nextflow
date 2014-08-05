/*
 * Copyright (c) 2013-2014, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2014, Paolo Di Tommaso and the respective authors.
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
import com.beust.jcommander.ParameterException
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class LauncherTest extends Specification {


    def testCommandLineOptions() {

        when:
        def opt = new Launcher('-daemon.x=1', '-daemon.y.z=2').parseMainArgs().options

        then:
        opt.daemonOptions.x == '1'
        opt.daemonOptions.'y.z'== '2'

    }


    def testVersion () {

        when:
        def launcher = new Launcher('-v').parseMainArgs()
        then:
        assert launcher.options.version

        when:
        launcher = new Launcher('-version').parseMainArgs()
        then:
        assert launcher.options.version
        assert launcher.fullVersion


    }

    def testHelp () {

        when:
        def launcher = new Launcher('-h').parseMainArgs()

        then:
        assert launcher.options.help

        when:
        launcher = new Launcher('help').parseMainArgs()
        then:
        launcher.command instanceof CmdHelp
        launcher.command.args == null

        when:
        launcher = new Launcher('help','xxx').parseMainArgs()
        then:
        launcher.command instanceof CmdHelp
        launcher.command.args == ['xxx']

    }

    def testInfo() {

        when:
        def launcher = new Launcher('info').parseMainArgs()
        then:
        launcher.command instanceof CmdInfo
        launcher.command.args == null

        when:
        launcher = new Launcher('info','xxx').parseMainArgs()
        then:
        launcher.command instanceof CmdInfo
        launcher.command.args == ['xxx']

    }

    def testPull() {

        when:
        def launcher = new Launcher('pull').parseMainArgs()
        then:
        thrown(ParameterException)

        when:
        launcher = new Launcher('pull','xxx').parseMainArgs()
        then:
        launcher.command instanceof CmdPull
        launcher.command.args == ['xxx']

    }


}
