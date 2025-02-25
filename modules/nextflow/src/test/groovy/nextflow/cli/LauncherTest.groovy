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

import java.nio.file.Files

import com.beust.jcommander.DynamicParameter
import com.beust.jcommander.Parameter
import spock.lang.Specification
import spock.lang.Unroll
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
        launcher.normalizeArgs('-x', '1', '-process.debug') == ['-x', '1', '-process.debug=true']
        launcher.normalizeArgs('-x', '1', '-process.debug', '-with-docker', 'ubuntu' ) == ['-x', '1', '-process.debug=true', '-with-docker','ubuntu']
        launcher.normalizeArgs('-x', '1', '-process.debug', '-123') == ['-x', '1', '-process.debug=-123' ]

        launcher.normalizeArgs('-x', '1', '-cluster.alpha','2', '3') == ['-x', '1', '-cluster.alpha=2', '3']
        launcher.normalizeArgs('-x', '1', '-cluster.debug') == ['-x', '1', '-cluster.debug=true']

        launcher.normalizeArgs('-x', '1', '-executor.alpha','2', '3') == ['-x', '1', '-executor.alpha=2', '3']
        launcher.normalizeArgs('-x', '1', '-executor.debug') == ['-x', '1', '-executor.debug=true']

        launcher.normalizeArgs('-x', '1', '-that.alpha','2', '3') == ['-x', '1', '-that.alpha','2', '3']

        launcher.normalizeArgs('run', 'file-name', '-a', '-b') == ['run','file-name', '-a', '-b']
        launcher.normalizeArgs('run', '-', '-a', '-b') == ['run','-stdin', '-a', '-b']
        launcher.normalizeArgs('run') == ['run']

        launcher.normalizeArgs('run','-with-cloudcache') == ['run', '-with-cloudcache', '-']
        launcher.normalizeArgs('run','-with-cloudcache', '-x') == ['run', '-with-cloudcache', '-', '-x']
        launcher.normalizeArgs('run','-with-cloudcache', 's3://foo/bar') == ['run', '-with-cloudcache','s3://foo/bar']
        
        launcher.normalizeArgs('run','-with-tower') == ['run', '-with-tower', '-']
        launcher.normalizeArgs('run','-with-tower', '-x') == ['run', '-with-tower', '-', '-x']
        launcher.normalizeArgs('run','-with-tower', 'foo.com') == ['run', '-with-tower','foo.com']

        launcher.normalizeArgs('run','-with-wave') == ['run', '-with-wave', '-']
        launcher.normalizeArgs('run','-with-wave', '-x') == ['run', '-with-wave', '-', '-x']
        launcher.normalizeArgs('run','-with-wave', 'foo.com') == ['run', '-with-wave','foo.com']

        launcher.normalizeArgs('run','-with-trace') == ['run', '-with-trace','-']
        launcher.normalizeArgs('run','-with-trace', '-x') == ['run', '-with-trace','-', '-x']
        launcher.normalizeArgs('run','-with-trace', 'file.x') == ['run', '-with-trace','file.x']

        launcher.normalizeArgs('run','-with-report') == ['run', '-with-report','-']
        launcher.normalizeArgs('run','-with-report', '-x') == ['run', '-with-report','-', '-x']
        launcher.normalizeArgs('run','-with-report', 'file.x') == ['run', '-with-report','file.x']

        launcher.normalizeArgs('run','-with-timeline') == ['run', '-with-timeline','-']
        launcher.normalizeArgs('run','-with-timeline', '-x') == ['run', '-with-timeline','-', '-x']
        launcher.normalizeArgs('run','-with-timeline', 'file.x') == ['run', '-with-timeline','file.x']

        launcher.normalizeArgs('run','-with-dag') == ['run', '-with-dag','-']
        launcher.normalizeArgs('run','-with-dag', '-x') == ['run', '-with-dag','-', '-x']
        launcher.normalizeArgs('run','-with-dag', 'file.dot') == ['run', '-with-dag','file.dot']

        launcher.normalizeArgs('run','-with-docker') == ['run', '-with-docker','-']
        launcher.normalizeArgs('run','-with-docker', '-x') == ['run', '-with-docker','-', '-x']
        launcher.normalizeArgs('run','-with-docker', 'busybox') == ['run', '-with-docker','busybox']

        launcher.normalizeArgs('run','-with-podman') == ['run', '-with-podman','-']
        launcher.normalizeArgs('run','-with-podman', '-x') == ['run', '-with-podman','-', '-x']
        launcher.normalizeArgs('run','-with-podman', 'busybox') == ['run', '-with-podman','busybox']

        launcher.normalizeArgs('run','-with-singularity') == ['run', '-with-singularity','-']
        launcher.normalizeArgs('run','-with-singularity', '-x') == ['run', '-with-singularity','-', '-x']
        launcher.normalizeArgs('run','-with-singularity', 'busybox') == ['run', '-with-singularity','busybox']

        launcher.normalizeArgs('run','-with-charliecloud') == ['run', '-with-charliecloud','-']
        launcher.normalizeArgs('run','-with-charliecloud', '-x') == ['run', '-with-charliecloud','-', '-x']
        launcher.normalizeArgs('run','-with-charliecloud', 'busybox') == ['run', '-with-charliecloud','busybox']

        launcher.normalizeArgs('run','-with-conda') == ['run', '-with-conda','-']
        launcher.normalizeArgs('run','-with-conda', '-x') == ['run', '-with-conda','-', '-x']
        launcher.normalizeArgs('run','-with-conda', 'busybox') == ['run', '-with-conda','busybox']

        launcher.normalizeArgs('run','-with-spack') == ['run', '-with-spack','-']
        launcher.normalizeArgs('run','-with-spack', '-x') == ['run', '-with-spack','-', '-x']
        launcher.normalizeArgs('run','-with-spack', 'busybox') == ['run', '-with-spack','busybox']

        launcher.normalizeArgs('run','-dump-channels') == ['run', '-dump-channels','*']
        launcher.normalizeArgs('run','-dump-channels', '-x') == ['run', '-dump-channels','*', '-x']
        launcher.normalizeArgs('run','-dump-channels', 'foo,bar') == ['run', '-dump-channels','foo,bar']

        launcher.normalizeArgs('run','-with-notification', 'paolo@yo.com') == ['run', '-with-notification','paolo@yo.com']
        launcher.normalizeArgs('run','-with-notification') == ['run', '-with-notification','true']
        launcher.normalizeArgs('run','-with-notification', '-x') == ['run', '-with-notification','true', '-x']

        launcher.normalizeArgs('run','-with-fusion', 'false') == ['run', '-with-fusion','false']
        launcher.normalizeArgs('run','-with-fusion') == ['run', '-with-fusion','true']
        launcher.normalizeArgs('run','-with-fusion', '-x') == ['run', '-with-fusion','true', '-x']

        launcher.normalizeArgs('run','-N', 'paolo@yo.com') == ['run', '-N','paolo@yo.com']
        launcher.normalizeArgs('run','-N') == ['run', '-N','true']
        launcher.normalizeArgs('run','-N', '-x') == ['run', '-N','true', '-x']

        launcher.normalizeArgs('run','-syslog', 'host.com') == ['run', '-syslog','host.com']
        launcher.normalizeArgs('run','-syslog') == ['run', '-syslog','localhost']
        launcher.normalizeArgs('run','-syslog', '-x') == ['run', '-syslog','localhost', '-x']

        launcher.normalizeArgs('run','-ansi-log', '-x') == ['run', '-ansi-log','true', '-x']
        launcher.normalizeArgs('run','-ansi-log', 'true', '-x') == ['run', '-ansi-log','true', '-x']
        launcher.normalizeArgs('run','-ansi-log', 'false', '-x') == ['run', '-ansi-log','false', '-x']

        launcher.normalizeArgs('run','-stub', '-x') == ['run', '-stub','true', '-x']
        launcher.normalizeArgs('run','-stub-run', '-x') == ['run', '-stub-run','true', '-x']
        launcher.normalizeArgs('run','-stub-run', 'true', '-x') == ['run', '-stub-run', 'true', '-x']
        launcher.normalizeArgs('run','-stub-run', 'false', '-x') == ['run', '-stub-run', 'false', '-x']

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

        when:
        Launcher.setProxy('FTP', [FTP_PROXY: 'delta.com:7566'])
        then:
        System.getProperty('ftp.proxyHost') == 'delta.com'
        System.getProperty('ftp.proxyPort') == '7566'

        when:
        Launcher.setProxy('ftp', [ftp_proxy: 'epsilon.com:6658'])
        then:
        System.getProperty('ftp.proxyHost') == 'epsilon.com'
        System.getProperty('ftp.proxyPort') == '6658'
    }

    @RestoreSystemProperties
    def 'should setup proxy properties and configure the network authenticator'() {

        when:
        Launcher.setProxy('HTTP', [HTTP_PROXY: 'http://alphauser:alphapass@alpha.com:333'])
        PasswordAuthentication auth = Authenticator.requestPasswordAuthentication(
            'alpha.com', null, 333, 'http', null, null, null, Authenticator.RequestorType.PROXY
        )
        then:
        System.getProperty('http.proxyHost') == 'alpha.com'
        System.getProperty('http.proxyPort') == '333'
        and:
        auth.getUserName() == 'alphauser'
        auth.getPassword() == 'alphapass'.toCharArray()

        when:
        Launcher.setProxy('http', [http_proxy: 'http://gammauser:gammapass@gamma.com:444'])
        auth = Authenticator.requestPasswordAuthentication(
            'gamma.com', null, 444, 'http', null, null, null, Authenticator.RequestorType.PROXY
        )
        then:
        System.getProperty('http.proxyHost') == 'gamma.com'
        System.getProperty('http.proxyPort') == '444'
        and:
        auth.getUserName() == 'gammauser'
        auth.getPassword() == 'gammapass'.toCharArray()

        when:
        Launcher.setProxy('HTTPS', [HTTPS_PROXY: 'https://betauser:betapass@beta.com:5466'])
        auth = Authenticator.requestPasswordAuthentication(
            'beta.com', null, 5466, 'HTTPS', null, null, null, Authenticator.RequestorType.PROXY
        )
        then:
        System.getProperty('https.proxyHost') == 'beta.com'
        System.getProperty('https.proxyPort') == '5466'
        and:
        auth.getUserName() == 'betauser'
        auth.getPassword() == 'betapass'.toCharArray()

        when:
        Launcher.setProxy('https', [https_proxy: 'https://zetauser:zetapass@zeta.com:6646'])
        auth = Authenticator.requestPasswordAuthentication(
            'zeta.com', null, 6646, 'https', null, null, null, Authenticator.RequestorType.PROXY
        )
        then:
        System.getProperty('https.proxyHost') == 'zeta.com'
        System.getProperty('https.proxyPort') == '6646'
        and:
        auth.getUserName() == 'zetauser'
        auth.getPassword() == 'zetapass'.toCharArray()

        when:
        Launcher.setProxy('FTP', [FTP_PROXY: 'ftp://deltauser:deltapass@delta.com:7566'])
        auth = Authenticator.requestPasswordAuthentication(
            'delta.com', null, 7566, 'ftp', null, null, null, Authenticator.RequestorType.PROXY
        )
        then:
        System.getProperty('ftp.proxyHost') == 'delta.com'
        System.getProperty('ftp.proxyPort') == '7566'
        and:
        auth.getUserName() == 'deltauser'
        auth.getPassword() == 'deltapass'.toCharArray()

        when:
        Launcher.setProxy('ftp', [ftp_proxy: 'ftp://epsilonuser:epsilonpass@epsilon.com:6658'])
        auth = Authenticator.requestPasswordAuthentication(
            'epsilon.com', null, 6658, 'ftp', null, null, null, Authenticator.RequestorType.PROXY
        )
        then:
        System.getProperty('ftp.proxyHost') == 'epsilon.com'
        System.getProperty('ftp.proxyPort') == '6658'
        and:
        auth.getUserName() == 'epsilonuser'
        auth.getPassword() == 'epsilonpass'.toCharArray()
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

    @RestoreSystemProperties
    @Unroll
    def 'should set http client timeout' () {
        when:
        Launcher.setHttpClientProperties(ENV)
        then:
        System.getProperty('jdk.httpclient.keepalive.timeout') == TIMEOUT
        and:
        System.getProperty('jdk.httpclient.connectionPoolSize') == POOLSIZE

        where:
        ENV                                             | TIMEOUT   | POOLSIZE
        [:]                                             | '10'      | null
        and:
        [NXF_JDK_HTTPCLIENT_KEEPALIVE_TIMEOUT: '1']     | '1'       | null
        [NXF_JDK_HTTPCLIENT_KEEPALIVE_TIMEOUT: '100']   | '100'     | null
        and:
        [NXF_JDK_HTTPCLIENT_CONNECTIONPOOLSIZE: '0']    | '10'      | '0'
        [NXF_JDK_HTTPCLIENT_CONNECTIONPOOLSIZE: '99']   | '10'      | '99'
    }

    def 'should make cli' () {
        given:
        def launcher = new Launcher()
        expect:
        launcher.makeCli('nextflow','run','foo.nf') == 'nextflow run foo.nf'
        launcher.makeCli('nextflow','run','foo.nf', '*.txt') == "nextflow run foo.nf '*.txt'"
        launcher.makeCli('/this/that/nextflow run foo.nf','run','foo.nf', 'a{1,2}.z') == "nextflow run foo.nf 'a{1,2}.z'"
        launcher.makeCli('/this/that/launch run bar.nf','run','bar.nf') == '/this/that/launch run bar.nf'
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
