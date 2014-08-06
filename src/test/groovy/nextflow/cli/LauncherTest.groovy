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
        def launcher = new Launcher('pull','xxx').parseMainArgs()
        then:
        launcher.command instanceof CmdPull
        launcher.command.args == ['xxx']

    }


    def testNormalizeCmdline () {

        expect:
        Launcher.normalizeArgs('a','-bb','-ccc','dddd') == ['a','-bb','-ccc','dddd']
        Launcher.normalizeArgs('a','-bb','-ccc','-resume', 'last') == ['a','-bb','-ccc','-resume','last']
        Launcher.normalizeArgs('a','-bb','-ccc','-resume') == ['a','-bb','-ccc','-resume','last']
        Launcher.normalizeArgs('a','-bb','-ccc','-resume','1d2c942a-345d-420b-b7c7-18d90afc6c33', 'zzz') == ['a','-bb','-ccc','-resume','1d2c942a-345d-420b-b7c7-18d90afc6c33', 'zzz']

        Launcher.normalizeArgs('x','-test') == ['x','-test','%all']
        Launcher.normalizeArgs('x','-test','alpha') == ['x','-test','alpha']
        Launcher.normalizeArgs('x','-test','-other') == ['x','-test','%all','-other']

        Launcher.normalizeArgs('--alpha=1') == ['--alpha=1']
        Launcher.normalizeArgs('--alpha','1') == ['--alpha=1']
        Launcher.normalizeArgs('-x', '1', 'script.nf', '--long', 'v1', '--more', 'v2', '--flag') == ['-x','1','script.nf','--long=v1','--more=v2','--flag=true']

        Launcher.normalizeArgs('-x', '1', '-process.alpha','2', '3') == ['-x', '1', '-process.alpha=2', '3']
        Launcher.normalizeArgs('-x', '1', '-process.echo') == ['-x', '1', '-process.echo=true']

        Launcher.normalizeArgs('-x', '1', '-daemon.alpha','2', '3') == ['-x', '1', '-daemon.alpha=2', '3']
        Launcher.normalizeArgs('-x', '1', '-daemon.echo') == ['-x', '1', '-daemon.echo=true']

        Launcher.normalizeArgs('-x', '1', '-executor.alpha','2', '3') == ['-x', '1', '-executor.alpha=2', '3']
        Launcher.normalizeArgs('-x', '1', '-executor.echo') == ['-x', '1', '-executor.echo=true']

        Launcher.normalizeArgs('-x', '1', '-that.alpha','2', '3') == ['-x', '1', '-that.alpha','2', '3']

        Launcher.normalizeArgs('run', 'file-name', '-a', '-b') == ['run','file-name', '-a', '-b']
        Launcher.normalizeArgs('run', '-', '-a', '-b') == ['run','-stdin', '-a', '-b']
        Launcher.normalizeArgs('run') == ['run']

    }



}
