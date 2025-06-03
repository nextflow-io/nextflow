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

package nextflow.pixi

import java.nio.file.Files
import java.nio.file.Paths

import nextflow.util.Duration
import spock.lang.IgnoreIf
import spock.lang.PendingFeature
import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Edmund Miller <edmund.miller@seqera.io>
 */
class PixiConfigTest extends Specification {

    @Unroll
    def 'should check enabled flag'() {
        given:
        def pixi = new PixiConfig(CONFIG, ENV)
        expect:
        pixi.isEnabled() == EXPECTED

        where:
        EXPECTED    | CONFIG            | ENV
        false       | [:]               | [:]
        false       | [enabled: false]  | [:]
        true        | [enabled: true]   | [:]
        and:
        false       | [:]               | [NXF_PIXI_ENABLED: 'false']
        true        | [:]               | [NXF_PIXI_ENABLED: 'true']
        false       | [enabled: false]  | [NXF_PIXI_ENABLED: 'true']  // <-- config has priority
        true        | [enabled: true]   | [NXF_PIXI_ENABLED: 'true']
    }

    def 'should return create timeout'() {
        given:
        def CONFIG = [createTimeout: '30min']
        def pixi = new PixiConfig(CONFIG, [:])

        expect:
        pixi.createTimeout() == Duration.of('30min')
    }

    def 'should return null create timeout when not specified'() {
        given:
        def pixi = new PixiConfig([:], [:])

        expect:
        pixi.createTimeout() == null
    }

    def 'should return create options'() {
        given:
        def CONFIG = [createOptions: '--verbose --no-lock-update']
        def pixi = new PixiConfig(CONFIG, [:])

        expect:
        pixi.createOptions() == '--verbose --no-lock-update'
    }

    def 'should return null create options when not specified'() {
        given:
        def pixi = new PixiConfig([:], [:])

        expect:
        pixi.createOptions() == null
    }

    def 'should return cache directory'() {
        given:
        def CONFIG = [cacheDir: '/my/cache/dir']
        def pixi = new PixiConfig(CONFIG, [:])

        expect:
        pixi.cacheDir() == Paths.get('/my/cache/dir')
    }

    def 'should return null cache directory when not specified'() {
        given:
        def pixi = new PixiConfig([:], [:])

        expect:
        pixi.cacheDir() == null
    }

    def 'should handle boolean values for enabled flag'() {
        given:
        def pixi = new PixiConfig([enabled: true], [:])

        expect:
        pixi.isEnabled() == true

        when:
        pixi = new PixiConfig([enabled: false], [:])

        then:
        pixi.isEnabled() == false
    }

    def 'should handle string values for enabled flag'() {
        given:
        def pixi = new PixiConfig([enabled: 'true'], [:])

        expect:
        pixi.isEnabled() == true

        when:
        pixi = new PixiConfig([enabled: 'false'], [:])

        then:
        pixi.isEnabled() == false
    }

    def 'should inherit from LinkedHashMap'() {
        given:
        def CONFIG = [enabled: true, createTimeout: '10min', customOption: 'value']
        def pixi = new PixiConfig(CONFIG, [:])

        expect:
        pixi instanceof LinkedHashMap
        pixi.enabled == true
        pixi.createTimeout == '10min'
        pixi.customOption == 'value'
        pixi.size() == 3
    }

    // ==== Integration Tests (require actual Pixi installation) ====

    /**
     * Check if Pixi is installed and available on the system
     */
    private static boolean hasPixiInstalled() {
        try {
            def process = new ProcessBuilder('pixi', '--version').start()
            process.waitFor()
            return process.exitValue() == 0
        } catch (Exception e) {
            return false
        }
    }

    @IgnoreIf({ !hasPixiInstalled() })
    def 'should verify pixi version command works'() {
        when:
        def process = new ProcessBuilder('pixi', '--version').start()
        def exitCode = process.waitFor()
        def output = process.inputStream.text.trim()

        then:
        exitCode == 0
        output.contains('pixi')
        output.matches(/pixi \d+\.\d+\.\d+.*/)
    }

    @IgnoreIf({ !hasPixiInstalled() })
    def 'should verify pixi info command works'() {
        when:
        def process = new ProcessBuilder('pixi', 'info').start()
        def exitCode = process.waitFor()
        def output = process.inputStream.text.trim()

        then:
        exitCode == 0
        output.contains('pixi')
        // Should contain platform information
        output.toLowerCase().contains('platform')
    }

    @IgnoreIf({ !hasPixiInstalled() })
    def 'should create and remove a temporary pixi environment'() {
        given:
        def tempDir = Files.createTempDirectory('pixi-test')
        def pixiToml = tempDir.resolve('pixi.toml')

        // Create a minimal pixi.toml
        pixiToml.text = '''
[project]
name = "test-env"
version = "0.1.0"
description = "Test environment"
channels = ["conda-forge"]
platforms = ["linux-64", "osx-64", "osx-arm64", "win-64"]

[dependencies]
python = ">=3.8"
'''.stripIndent()

        when:
        // Create the environment
        def createProcess = new ProcessBuilder('pixi', 'install')
            .directory(tempDir.toFile())
            .start()
        def createExitCode = createProcess.waitFor()

        then:
        createExitCode == 0
        tempDir.resolve('.pixi').exists()

        cleanup:
        tempDir?.deleteDir()
    }

    @IgnoreIf({ !hasPixiInstalled() })
    def 'should validate pixi environment creation with specific dependencies'() {
        given:
        def tempDir = Files.createTempDirectory('pixi-test-deps')
        def pixiToml = tempDir.resolve('pixi.toml')

        // Create a pixi.toml with specific dependencies commonly used in bioinformatics
        pixiToml.text = '''
[project]
name = "bio-test-env"
version = "0.1.0"
description = "Bioinformatics test environment"
channels = ["conda-forge", "bioconda"]
platforms = ["linux-64", "osx-64", "osx-arm64"]

[dependencies]
python = ">=3.9"
numpy = "*"
'''.stripIndent()

        when:
        // Install the environment
        def installProcess = new ProcessBuilder('pixi', 'install')
            .directory(tempDir.toFile())
            .start()
        def installExitCode = installProcess.waitFor()

        and:
        // List the installed packages
        def listProcess = new ProcessBuilder('pixi', 'list')
            .directory(tempDir.toFile())
            .start()
        def listExitCode = listProcess.waitFor()
        def listOutput = listProcess.inputStream.text

        then:
        installExitCode == 0
        listExitCode == 0
        tempDir.resolve('.pixi').exists()
        listOutput.contains('python')
        listOutput.contains('numpy')

        cleanup:
        tempDir?.deleteDir()
    }

    @IgnoreIf({ !hasPixiInstalled() })
    def 'should handle pixi environment with lock file'() {
        given:
        def tempDir = Files.createTempDirectory('pixi-test-lock')
        def pixiToml = tempDir.resolve('pixi.toml')

        pixiToml.text = '''
[project]
name = "lock-test-env"
version = "0.1.0"
description = "Lock file test environment"
channels = ["conda-forge"]
platforms = ["linux-64", "osx-64", "osx-arm64"]

[dependencies]
python = "3.11.*"
'''.stripIndent()

        when:
        // Create lock file
        def lockProcess = new ProcessBuilder('pixi', 'install')
            .directory(tempDir.toFile())
            .start()
        def lockExitCode = lockProcess.waitFor()

        then:
        lockExitCode == 0
        tempDir.resolve('pixi.lock').exists()
        tempDir.resolve('.pixi').exists()

        when:
        // Clean and reinstall from lock file
        tempDir.resolve('.pixi').deleteDir()
        def reinstallProcess = new ProcessBuilder('pixi', 'install')
            .directory(tempDir.toFile())
            .start()
        def reinstallExitCode = reinstallProcess.waitFor()

        then:
        reinstallExitCode == 0
        tempDir.resolve('.pixi').exists()

        cleanup:
        tempDir?.deleteDir()
    }

    @IgnoreIf({ !hasPixiInstalled() })
    def 'should run commands in pixi environment'() {
        given:
        def tempDir = Files.createTempDirectory('pixi-test-run')
        def pixiToml = tempDir.resolve('pixi.toml')

        pixiToml.text = '''
[project]
name = "run-test-env"
version = "0.1.0"
description = "Run command test environment"
channels = ["conda-forge"]
platforms = ["linux-64", "osx-64", "osx-arm64"]

[dependencies]
python = ">=3.9"

[tasks]
hello = "python -c 'print(\\\"Hello from Pixi!\\\")'"
'''.stripIndent()

        when:
        // Install the environment
        def installProcess = new ProcessBuilder('pixi', 'install')
            .directory(tempDir.toFile())
            .redirectErrorStream(true)
            .start()
        def installOutput = installProcess.inputStream.text
        def installExitCode = installProcess.waitFor()

        and:
        // Run the hello task
        def runProcess = new ProcessBuilder('pixi', 'run', 'hello')
            .directory(tempDir.toFile())
            .redirectErrorStream(true)
            .start()
        def runOutput = runProcess.inputStream.text
        def runExitCode = runProcess.waitFor()

        then:
        if (installExitCode != 0) {
            println "Pixi install failed with output:\n${installOutput}"
        }
        installExitCode == 0
        runExitCode == 0
        assert runOutput.contains('Hello from Pixi!') : "Pixi output was:\n${runOutput}"

        cleanup:
        tempDir?.deleteDir()
    }

    @IgnoreIf({ !hasPixiInstalled() })
    def 'should validate pixi global commands'() {
        when:
        // Test pixi global list (should work even if no global packages are installed)
        def globalListProcess = new ProcessBuilder('pixi', 'global', 'list').start()
        def globalListExitCode = globalListProcess.waitFor()

        then:
        globalListExitCode == 0

        when:
        // Test pixi search for a common package
        def searchProcess = new ProcessBuilder('pixi', 'search', 'python').start()
        def searchExitCode = searchProcess.waitFor()
        def searchOutput = searchProcess.inputStream.text

        then:
        searchExitCode == 0
        searchOutput.contains('python')
    }
}
