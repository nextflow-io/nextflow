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

import spock.lang.Specification

import java.nio.file.Files

import com.beust.jcommander.DynamicParameter
import com.beust.jcommander.Parameter
import spock.util.environment.RestoreSystemProperties
import test.OutputCapture
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class LauncherTest extends Specification {

    @org.junit.Rule
    OutputCapture capture = new OutputCapture()

    def 'should return `version` option' () {

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

    def 'should return `help` command' () {

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

    def 'should return `info` command'() {

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

    def 'should return `pull` command'() {

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

    def 'should return `clone` command'() {
        when:
        def launcher = new Launcher().parseMainArgs('clone','xxx', '-hub', 'bitbucket', '-user','xx:yy')
        then:
        launcher.command instanceof CmdClone
        launcher.command.args == ['xxx']
        launcher.command.hubProvider == 'bitbucket'
        launcher.command.hubUser == 'xx'
        launcher.command.hubPassword == 'yy'
    }


    def 'should return `run` command'() {
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
        launcher = new Launcher().parseMainArgs('run', 'script.nf', '--alpha', '0', '--omega', '9')
        then:
        launcher.command instanceof CmdRun
        launcher.command.params.'alpha' == '0'
        launcher.command.params.'omega' == '9'

    }


    def 'should normalise command line options' () {

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
        launcher.normalizeArgs('-x', '1', '-process.echo', '-with-docker', 'ubuntu' ) == ['-x', '1', '-process.echo=true', '-with-docker','ubuntu']
        launcher.normalizeArgs('-x', '1', '-process.echo', '-123') == ['-x', '1', '-process.echo=-123' ]

        launcher.normalizeArgs('-x', '1', '-cluster.alpha','2', '3') == ['-x', '1', '-cluster.alpha=2', '3']
        launcher.normalizeArgs('-x', '1', '-cluster.echo') == ['-x', '1', '-cluster.echo=true']

        launcher.normalizeArgs('-x', '1', '-executor.alpha','2', '3') == ['-x', '1', '-executor.alpha=2', '3']
        launcher.normalizeArgs('-x', '1', '-executor.echo') == ['-x', '1', '-executor.echo=true']

        launcher.normalizeArgs('-x', '1', '-that.alpha','2', '3') == ['-x', '1', '-that.alpha','2', '3']

        launcher.normalizeArgs('run', 'file-name', '-a', '-b') == ['run','file-name', '-a', '-b']
        launcher.normalizeArgs('run', '-', '-a', '-b') == ['run','-stdin', '-a', '-b']
        launcher.normalizeArgs('run') == ['run']

        launcher.normalizeArgs('run','-with-trace') == ['run', '-with-trace','trace.txt']
        launcher.normalizeArgs('run','-with-trace', '-x') == ['run', '-with-trace','trace.txt', '-x']
        launcher.normalizeArgs('run','-with-trace', 'file.x') == ['run', '-with-trace','file.x']

        launcher.normalizeArgs('run','-with-report') == ['run', '-with-report','report.html']
        launcher.normalizeArgs('run','-with-report', '-x') == ['run', '-with-report','report.html', '-x']
        launcher.normalizeArgs('run','-with-report', 'file.x') == ['run', '-with-report','file.x']

        launcher.normalizeArgs('run','-with-timeline') == ['run', '-with-timeline','timeline.html']
        launcher.normalizeArgs('run','-with-timeline', '-x') == ['run', '-with-timeline','timeline.html', '-x']
        launcher.normalizeArgs('run','-with-timeline', 'file.x') == ['run', '-with-timeline','file.x']

        launcher.normalizeArgs('run','-with-dag') == ['run', '-with-dag','dag.dot']
        launcher.normalizeArgs('run','-with-dag', '-x') == ['run', '-with-dag','dag.dot', '-x']
        launcher.normalizeArgs('run','-with-dag', 'file.dot') == ['run', '-with-dag','file.dot']

        launcher.normalizeArgs('run','-with-docker') == ['run', '-with-docker','-']
        launcher.normalizeArgs('run','-with-docker', '-x') == ['run', '-with-docker','-', '-x']
        launcher.normalizeArgs('run','-with-docker', 'busybox') == ['run', '-with-docker','busybox']

        launcher.normalizeArgs('run','-with-singularity') == ['run', '-with-singularity','-']
        launcher.normalizeArgs('run','-with-singularity', '-x') == ['run', '-with-singularity','-', '-x']
        launcher.normalizeArgs('run','-with-singularity', 'busybox') == ['run', '-with-singularity','busybox']

        launcher.normalizeArgs('run','-dump-channels') == ['run', '-dump-channels','*']
        launcher.normalizeArgs('run','-dump-channels', '-x') == ['run', '-dump-channels','*', '-x']
        launcher.normalizeArgs('run','-dump-channels', 'foo,bar') == ['run', '-dump-channels','foo,bar']

        launcher.normalizeArgs('run','-with-notification', 'paolo@yo.com') == ['run', '-with-notification','paolo@yo.com']
        launcher.normalizeArgs('run','-with-notification') == ['run', '-with-notification','true']
        launcher.normalizeArgs('run','-with-notification', '-x') == ['run', '-with-notification','true', '-x']

        launcher.normalizeArgs('run','-N', 'paolo@yo.com') == ['run', '-N','paolo@yo.com']
        launcher.normalizeArgs('run','-N') == ['run', '-N','true']
        launcher.normalizeArgs('run','-N', '-x') == ['run', '-N','true', '-x']

        launcher.normalizeArgs('run','-K', 'true') == ['run', '-K','true']
        launcher.normalizeArgs('run','-K') == ['run', '-K','true']
        launcher.normalizeArgs('run','-K', '-x') == ['run', '-K','true', '-x']

        launcher.normalizeArgs('run','-with-k8s', 'true') == ['run', '-with-k8s','true']
        launcher.normalizeArgs('run','-with-k8s') == ['run', '-with-k8s','true']
        launcher.normalizeArgs('run','-with-k8s', '-x') == ['run', '-with-k8s','true', '-x']

        launcher.normalizeArgs('run','-syslog', 'host.com') == ['run', '-syslog','host.com']
        launcher.normalizeArgs('run','-syslog') == ['run', '-syslog','localhost']
        launcher.normalizeArgs('run','-syslog', '-x') == ['run', '-syslog','localhost', '-x']

        launcher.normalizeArgs( script.toAbsolutePath().toString(), '--x=1' ) == ['run', script.toAbsolutePath().toString(), '--x=1']

        launcher.normalizeArgs('--foo', '--bar x') == ['--foo=--bar x']

        cleanup:
        script?.delete()
    }

    def 'should validate isValue' () {
        expect:
        Launcher.isValue(STR) == EXPECTED
        
        where:
        STR                 | EXPECTED
        'foo'               | true
        '10'                | true
        '-10'               | true
        '-foo'              | false
        '--bar'             | false
        'x y'               | true
        '-x y'              | true
        '-x=y'              | false
    }


    def 'should parse proxy env variables'( ) {

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

    @RestoreSystemProperties
    def 'should setup proxy properties'() {

        when:
        Launcher.setProxy('HTTP', [HTTP_PROXY: 'alpha.com:333'])
        then:
        System.getProperty('http.proxyHost') == 'alpha.com'
        System.getProperty('http.proxyPort') == '333'

        when:
        Launcher.setProxy('http', [http_proxy: 'gamma.com:444'])
        then:
        System.getProperty('http.proxyHost') == 'gamma.com'
        System.getProperty('http.proxyPort') == '444'

        when:
        Launcher.setProxy('HTTPS', [HTTPS_PROXY: 'beta.com:5466'])
        then:
        System.getProperty('https.proxyHost') == 'beta.com'
        System.getProperty('https.proxyPort') == '5466'

        when:
        Launcher.setProxy('https', [https_proxy: 'zeta.com:6646'])
        then:
        System.getProperty('https.proxyHost') == 'zeta.com'
        System.getProperty('https.proxyPort') == '6646'

    }

    @RestoreSystemProperties
    def 'should set no proxy property' () {

        given:
        System.properties.remove('http.nonProxyHosts')
        
        when:
        Launcher.setNoProxy(ENV)
        then:
        System.getProperty('http.nonProxyHosts') == EXPECTED

        where:
        ENV                         | EXPECTED
        [:]                         | null
        [no_proxy: 'localhost' ]    | 'localhost'
        [NO_PROXY: '127.0.0.1' ]    | '127.0.0.1'
        [NO_PROXY:'localhost,127.0.0.1,.localdomain.com']  | 'localhost|127.0.0.1|.localdomain.com'

    }

    def 'should print Parameter and DynamicParameter annotation'() {

        given:
        def launcher = new Launcher()
        when:
        launcher.printOptions(Opts)
        then:
        capture.toString() == '''
            Options:
              -D
                 Set JVM properties
              -c, -config
                 Add the specified file to configuration set
              -log
                 Set nextflow log file
            '''
                .stripIndent().leftTrim()
    }

    def 'should print list of commands'() {
        given:
        def launcher = new Launcher()
        when:
        launcher.printCommands( [new CmdInfo(), new CmdRun(), new CmdList()] )
        then:
        capture.toString() == '''
                Commands:
                  info   Print project and system runtime information
                  list   List all downloaded projects
                  run    Execute a pipeline project

               '''
                .stripIndent()
    }

    static class Opts {

        @Parameter(names=['-log'], description = 'Set nextflow log file')
        String opt1

        @Parameter(names=['-c','-config'], description = 'Add the specified file to configuration set')
        String opt2

        @DynamicParameter(names = ['-D'], description = 'Set JVM properties' )
        Map opt3

        @Parameter(names=['hidden'], hidden = true)
        String opt4

    }
}
