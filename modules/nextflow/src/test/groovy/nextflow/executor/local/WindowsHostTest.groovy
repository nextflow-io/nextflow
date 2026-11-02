/*
 * Copyright 2013-2026, Seqera Labs
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

package nextflow.executor.local

import java.nio.file.Files
import java.nio.file.Path

import nextflow.Session
import nextflow.container.DockerConfig
import nextflow.executor.BashWrapperBuilder
import nextflow.executor.SimpleFileCopyStrategy
import nextflow.file.FileHelper
import nextflow.processor.TaskBean
import nextflow.secret.LocalSecretsProvider
import nextflow.secret.SecretImpl
import nextflow.util.Escape
import nextflow.util.PosixPaths
import nextflow.util.SysHelper
import spock.lang.Requires
import spock.lang.Specification

/**
 * Checks the support for running Nextflow natively on Windows using the actual
 * Windows file system paths. These tests only run on a Windows host, see the
 * `windows.yml` GitHub Actions workflow.
 */
@Requires({ os.windows })
class WindowsHostTest extends Specification {

    def setupSpec() {
        new Session()
    }

    def 'should detect the windows host' () {
        expect:
        SysHelper.isWindows()
        FileHelper.windowsHost
        BashWrapperBuilder.systemOsName.startsWith('Windows')
    }

    def 'should skip the signals not available on windows' () {
        expect:
        !Session.handleSignal('HUP', { sig -> } as sun.misc.SignalHandler)
    }

    def 'should store secrets on windows' () {
        given:
        def folder = Files.createTempDirectory('test')
        def storeFile = folder.resolve('secrets/store.json')
        Files.createDirectories(storeFile.parent)
        def provider = new LocalSecretsProvider(storeFile: storeFile)

        when:
        provider.storeSecrets([new SecretImpl('foo', 'bar')])
        then:
        storeFile.text.contains('"bar"')

        cleanup:
        folder?.deleteDir()
    }

    def 'should map windows paths to container paths' () {
        expect:
        PosixPaths.toPosix(Path.of('C:\\Users\\me\\work')) == '/c/Users/me/work'
        PosixPaths.toPosix(Path.of('D:/My Data/reads.fq')) == '/d/My Data/reads.fq'
        PosixPaths.toPosix(Path.of('C:\\')) == '/c'
    }

    def 'should parse drive paths with forward slashes as local paths' () {
        expect:
        FileHelper.asPath('C:/data/file.txt') == Path.of('C:\\data\\file.txt')
        FileHelper.asPath('C:\\data\\file.txt') == Path.of('C:\\data\\file.txt')
    }

    def 'should enable the docker wrapper mode' () {
        given:
        def executor = new LocalExecutor()
        executor.session = Mock(Session) {
            getConfig() >> [:]
            getContainerConfig() >> CONFIG
        }

        expect:
        executor.isDockerWrapperEnabled() == EXPECTED
        executor.isContainerNative() == EXPECTED

        where:
        CONFIG                              | EXPECTED
        new DockerConfig(enabled: true)     | true
        new DockerConfig(enabled: false)    | false
    }

    def 'should mount the task paths in the container' () {
        given:
        def bean = new TaskBean(
            name: 'foo',
            workDir: Path.of('C:\\work\\ab\\cdef'),
            containerImage: 'ubuntu:24.04',
            containerConfig: new DockerConfig(enabled: true),
            inputFiles: [
                'a.txt': Path.of('D:\\My Data\\a.txt'),
                'b.txt': Path.of('C:\\work\\xy\\z123\\b.txt') ],
            binDirs: [Path.of('C:\\project\\bin')] )

        expect:
        new DockerWrapperLauncher(bean, 'nxf-1234').runCommand() == [
            'docker', 'run', '--rm',
            '-v', 'C:\\project\\bin:/c/project/bin',
            '-v', 'C:\\work\\ab\\cdef:/c/work/ab/cdef',
            '-v', 'C:\\work\\xy\\z123:/c/work/xy/z123',
            '-v', 'D:\\My Data:/d/My Data',
            '-w', '/c/work/ab/cdef',
            '--name', 'nxf-1234',
            'ubuntu:24.04',
            '/bin/bash', '-ue', '.command.run' ]
    }

    def 'should stage input files using container paths' () {
        given:
        def workDir = Path.of('C:\\work\\ab\\cdef')
        def strategy = new SimpleFileCopyStrategy(new TaskBean(workDir: workDir, targetDir: workDir, stageInMode: MODE))

        expect:
        strategy.getStageInputFilesScript(['a.txt': Path.of(SOURCE)]) == "rm -f a.txt\n$EXPECTED"

        where:
        MODE        | SOURCE                        | EXPECTED
        'symlink'   | 'D:\\My Data\\a.txt'          | 'ln -sfn /d/My\\ Data/a.txt a.txt'
        'copy'      | 'D:\\data\\a.txt'             | 'cp -fRL /d/data/a.txt a.txt'
        'rellink'   | 'C:\\work\\xy\\z123\\a.txt'   | 'ln -sfn ../../xy/z123/a.txt a.txt'
    }

    def 'should create the task wrapper with container paths and unix line endings' () {
        given:
        def folder = Files.createTempDirectory('test')
        def bean = new TaskBean(
            name: 'foo',
            workDir: folder,
            script: 'echo Hello\r\necho world',
            containerImage: 'ubuntu:24.04',
            containerConfig: new DockerConfig(enabled: true),
            containerEnabled: true,
            containerNative: true,
            inputFiles: [:],
            outputFiles: [] )
        def builder = new BashWrapperBuilder(bean) {
            @Override
            protected String getSecretsEnv() { return null }
        }

        when:
        builder.build()
        def wrapper = folder.resolve('.command.run').text
        def script = folder.resolve('.command.sh').text
        then:
        !wrapper.contains('\r')
        !script.contains('\r')
        wrapper.contains("/bin/bash -ue ${Escape.path(PosixPaths.toPosix(folder.resolve('.command.sh')))}")
        !(wrapper =~ /[A-Za-z]:\\/)

        cleanup:
        folder?.deleteDir()
    }

}
