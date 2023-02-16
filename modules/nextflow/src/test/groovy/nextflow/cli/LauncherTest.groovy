/*
 * Copyright 2020-2022, Seqera Labs
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

import picocli.CommandLine
import picocli.CommandLine.HelpCommand
import spock.lang.Specification
import spock.util.environment.RestoreSystemProperties
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class LauncherTest extends Specification {

    def parseArgs(String... args) {
        def launcher = new Launcher()
        def cmd = new CommandLine(launcher)
        cmd.parseArgs(args)

        return launcher
    }

    def getCommand(Launcher launcher) {
        launcher.spec.commandLine().getParseResult().subcommand().commandSpec().commandLine().getCommand()
    }

    def getArgs(Launcher launcher) {
        launcher.spec.commandLine().getParseResult().originalArgs()
    }

    def 'should return `version` option' () {

        when:
        def launcher = parseArgs('-v')
        then:
        assert launcher.options.version

        when:
        launcher = parseArgs('-version')
        then:
        assert launcher.options.fullVersion

    }

    def 'should return `help` command' () {

        when:
        def launcher = parseArgs('-h')
        then:
        assert launcher.help

        when:
        launcher = parseArgs('help')
        then:
        getCommand(launcher) instanceof HelpCommand

        when:
        launcher = parseArgs('help', 'xxx')
        then:
        getCommand(launcher) instanceof HelpCommand

    }

    def 'should return `info` command'() {

        when:
        def launcher = parseArgs('info')
        def command = getCommand(launcher)
        then:
        command instanceof CmdInfo
        command.pipeline == null

        when:
        launcher = parseArgs('info', 'xxx')
        command = getCommand(launcher)
        then:
        command instanceof CmdInfo
        command.pipeline == 'xxx'

    }

    def 'should return `pull` command'() {

        when:
        def launcher = parseArgs('pull', 'alpha')
        def command = getCommand(launcher)
        then:
        command instanceof CmdPull
        command.pipeline == 'alpha'

        when:
        launcher = parseArgs('pull', 'xxx', '-hub', 'bitbucket', '-user', 'xx:11')
        command = getCommand(launcher)
        then:
        command instanceof CmdPull
        command.pipeline == 'xxx'
        command.hubProvider == 'bitbucket'
        command.hubUser == 'xx'
        command.hubPassword == '11'

    }

    def 'should return `clone` command'() {
        when:
        def launcher = parseArgs('clone', 'xxx', '-hub', 'bitbucket', '-user', 'xx:yy')
        def command = getCommand(launcher)
        then:
        command instanceof CmdClone
        command.pipeline == 'xxx'
        command.hubProvider == 'bitbucket'
        command.hubUser == 'xx'
        command.hubPassword == 'yy'
    }


    def 'should return `run` command'() {
        when:
        def launcher = parseArgs('run', 'xxx', '-hub', 'bitbucket', '-user', 'xx:yy')
        def command = getCommand(launcher)
        then:
        command instanceof CmdRun
        command.pipeline == 'xxx'
        command.hubProvider == 'bitbucket'
        command.hubUser == 'xx'
        command.hubPassword == 'yy'

        when:
        launcher = parseArgs('run', 'alpha', '-hub', 'github')
        command = getCommand(launcher)
        then:
        command instanceof CmdRun
        command.pipeline == 'alpha'
        command.hubProvider == 'github'

        when:
        launcher = parseArgs('run', 'script.nf', 'arg1', 'arg2', '--', '--alpha', '0', '--omega', '9')
        command = getCommand(launcher)
        then:
        command instanceof CmdRun
        command.args == ['arg1', 'arg2']
        command.params['alpha'] == '0'
        command.params['omega'] == '9'

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

    def 'should make cli' () {
        given:
        def launcher = new Launcher()
        expect:
        launcher.makeCli('nextflow', 'run', 'foo.nf') == 'nextflow run foo.nf'
        launcher.makeCli('nextflow', 'run', 'foo.nf', '*.txt') == "nextflow run foo.nf '*.txt'"
        launcher.makeCli('/this/that/nextflow run foo.nf', 'run', 'foo.nf', 'a{1,2}.z') == "nextflow run foo.nf 'a{1,2}.z'"
        launcher.makeCli('/this/that/launch run bar.nf', 'run', 'bar.nf') == '/this/that/launch run bar.nf'
    }

}
