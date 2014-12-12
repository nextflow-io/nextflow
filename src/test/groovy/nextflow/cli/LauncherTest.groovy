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

import java.nio.file.Files

import com.beust.jcommander.ParameterException
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class LauncherTest extends Specification {


    def testVersion () {

        when:
        def launcher = new Launcher().parseMainArgs('-v')
        then:
        assert launcher.options.version

        when:
        launcher = new Launcher().parseMainArgs('-version')
        then:
        assert launcher.options.version
        assert launcher.fullVersion


    }

    def testHelp () {

        when:
        def launcher = new Launcher().parseMainArgs('-h')

        then:
        assert launcher.options.help

        when:
        launcher = new Launcher().parseMainArgs('help')
        then:
        launcher.command instanceof CmdHelp
        launcher.command.args == null

        when:
        launcher = new Launcher().parseMainArgs('help','xxx')
        then:
        launcher.command instanceof CmdHelp
        launcher.command.args == ['xxx']

    }

    def testInfo() {

        when:
        def launcher = new Launcher().parseMainArgs('info')
        then:
        launcher.command instanceof CmdInfo
        launcher.command.args == null

        when:
        launcher = new Launcher().parseMainArgs('info','xxx')
        then:
        launcher.command instanceof CmdInfo
        launcher.command.args == ['xxx']

    }

    def testPull() {

        when:
        def launcher = new Launcher().parseMainArgs('pull','alpha')
        then:
        launcher.command instanceof CmdPull
        launcher.command.args == ['alpha']

        when:
        launcher = new Launcher().parseMainArgs('pull','xxx', '-hub', 'bitbucket', '-user','xx:11')
        then:
        launcher.command instanceof CmdPull
        launcher.command.args == ['xxx']
        launcher.command.hubProvider == 'bitbucket'
        launcher.command.hubUser == 'xx'
        launcher.command.hubPassword == '11'

    }

    def testClone() {
        when:
        def launcher = new Launcher().parseMainArgs('clone','xxx', '-hub', 'bitbucket', '-user','xx:yy')
        then:
        launcher.command instanceof CmdClone
        launcher.command.args == ['xxx']
        launcher.command.hubProvider == 'bitbucket'
        launcher.command.hubUser == 'xx'
        launcher.command.hubPassword == 'yy'
    }


    def testRun() {
        when:
        def launcher = new Launcher().parseMainArgs('run','xxx', '-hub', 'bitbucket', '-user','xx:yy')
        then:
        launcher.command instanceof CmdRun
        launcher.command.args == ['xxx']
        launcher.command.hubProvider == 'bitbucket'
        launcher.command.hubUser == 'xx'
        launcher.command.hubPassword == 'yy'

        when:
        launcher = new Launcher().parseMainArgs('run','alpha', '-hub', 'github')
        then:
        launcher.command instanceof CmdRun
        launcher.command.args == ['alpha']
        launcher.command.hubProvider == 'github'

        when:
        new Launcher().parseMainArgs('run','alpha', '-hub', 'xyz')
        then:
        thrown(ParameterException)
    }


    def testNormalizeCmdline () {

        given:
        def script = Files.createTempFile('file',null)
        def launcher = [:] as Launcher
        launcher.allCommands = [ new CmdRun(), new CmdInfo() ]

        expect:
        launcher.normalizeArgs('a','-bb','-ccc','dddd') == ['a','-bb','-ccc','dddd']
        launcher.normalizeArgs('a','-bb','-ccc','-resume', 'last') == ['a','-bb','-ccc','-resume','last']
        launcher.normalizeArgs('a','-bb','-ccc','-resume') == ['a','-bb','-ccc','-resume','last']
        launcher.normalizeArgs('a','-bb','-ccc','-resume','1d2c942a-345d-420b-b7c7-18d90afc6c33', 'zzz') == ['a','-bb','-ccc','-resume','1d2c942a-345d-420b-b7c7-18d90afc6c33', 'zzz']

        launcher.normalizeArgs('x','-test') == ['x','-test','%all']
        launcher.normalizeArgs('x','-test','alpha') == ['x','-test','alpha']
        launcher.normalizeArgs('x','-test','-other') == ['x','-test','%all','-other']

        launcher.normalizeArgs('--alpha=1') == ['--alpha=1']
        launcher.normalizeArgs('--alpha','1') == ['--alpha=1']
        launcher.normalizeArgs('run','--x') == ['run', '--x=true']
        launcher.normalizeArgs('run','--x','--y') == ['run', '--x=true', '--y=true']
        launcher.normalizeArgs('run','--x','--y', '-1', '--z') == ['run', '--x=true', '--y=-1', '--z=true']

        launcher.normalizeArgs('-x', '1', 'script.nf', '--long', 'v1', '--more', 'v2', '--flag') == ['-x','1','script.nf','--long=v1','--more=v2','--flag=true']

        launcher.normalizeArgs('-x', '1', '-process.alpha','2', '3') == ['-x', '1', '-process.alpha=2', '3']
        launcher.normalizeArgs('-x', '1', '-process.echo') == ['-x', '1', '-process.echo=true']

        launcher.normalizeArgs('-x', '1', '-cluster.alpha','2', '3') == ['-x', '1', '-cluster.alpha=2', '3']
        launcher.normalizeArgs('-x', '1', '-cluster.echo') == ['-x', '1', '-cluster.echo=true']

        launcher.normalizeArgs('-x', '1', '-executor.alpha','2', '3') == ['-x', '1', '-executor.alpha=2', '3']
        launcher.normalizeArgs('-x', '1', '-executor.echo') == ['-x', '1', '-executor.echo=true']

        launcher.normalizeArgs('-x', '1', '-that.alpha','2', '3') == ['-x', '1', '-that.alpha','2', '3']

        launcher.normalizeArgs('run', 'file-name', '-a', '-b') == ['run','file-name', '-a', '-b']
        launcher.normalizeArgs('run', '-', '-a', '-b') == ['run','-stdin', '-a', '-b']
        launcher.normalizeArgs('run') == ['run']

        launcher.normalizeArgs('run','-with-drmaa') == ['run', '-with-drmaa','-']
        launcher.normalizeArgs('run','-with-drmaa', '-x') == ['run', '-with-drmaa','-', '-x']
        launcher.normalizeArgs('run','-with-drmaa', 'X') == ['run', '-with-drmaa','X']

        launcher.normalizeArgs('run','-with-trace') == ['run', '-with-trace','trace.csv']
        launcher.normalizeArgs('run','-with-trace', '-x') == ['run', '-with-trace','trace.csv', '-x']
        launcher.normalizeArgs('run','-with-trace', 'file.x') == ['run', '-with-trace','file.x']

        launcher.normalizeArgs('run','-with-docker') == ['run', '-with-docker','-']
        launcher.normalizeArgs('run','-with-docker', '-x') == ['run', '-with-docker','-', '-x']
        launcher.normalizeArgs('run','-with-docker', 'busybox') == ['run', '-with-docker','busybox']

        launcher.normalizeArgs( script.toAbsolutePath().toString(), '--x=1' ) == ['run', script.toAbsolutePath().toString(), '--x=1']


        cleanup:
        script?.delete()
    }


    def testParseProxy( ) {

        expect:
        Launcher.parseProxy(null) == []
        Launcher.parseProxy('http://domain') == ['domain']
        Launcher.parseProxy('http://domain:333') == ['domain', '333']
        Launcher.parseProxy('http://10.20.30.40') == ['10.20.30.40']
        Launcher.parseProxy('http://10.20.30.40:333') == ['10.20.30.40', '333']
        Launcher.parseProxy('http://10.20.30.40:333/some/path') == ['10.20.30.40', '333']

        Launcher.parseProxy('foo') == ['foo']
        Launcher.parseProxy('foo:123') == ['foo','123']

    }

    def testSetupProxy() {

        when:
        Launcher.setProxy('http', [HTTP_PROXY: 'alpha.com:333'])
        then:
        System.getProperty('http.proxyHost') == 'alpha.com'
        System.getProperty('http.proxyPort') == '333'

        when:
        Launcher.setProxy('https', [HTTPS_PROXY: 'beta.com:5466'])
        then:
        System.getProperty('https.proxyHost') == 'beta.com'
        System.getProperty('https.proxyPort') == '5466'

    }


}
