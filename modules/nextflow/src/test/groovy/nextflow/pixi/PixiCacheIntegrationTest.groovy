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
import java.nio.file.Path
import java.nio.file.Paths

import nextflow.util.Duration
import spock.lang.IgnoreIf
import spock.lang.Specification

/**
 * Integration tests for PixiCache that require actual Pixi installation
 *
 * @author Edmund Miller <edmund.miller@seqera.io>
 */
class PixiCacheIntegrationTest extends Specification {

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
    def 'should create pixi environment from package specification'() {
        given:
        def tempDir = Files.createTempDirectory('pixi-cache-test')
        def config = new PixiConfig([cacheDir: tempDir.toString()], [:])
        def cache = new PixiCache(config)
        def ENV = 'python>=3.9'

        when:
        def prefix = cache.pixiPrefixPath(ENV)

        then:
        prefix != null
        prefix.parent == tempDir
        prefix.fileName.toString().startsWith('env-')

        cleanup:
        tempDir?.deleteDir()
    }

    @IgnoreIf({ !hasPixiInstalled() })
    def 'should create pixi environment from TOML file'() {
        given:
        def tempDir = Files.createTempDirectory('pixi-cache-toml-test')
        def config = new PixiConfig([cacheDir: tempDir.toString()], [:])
        def cache = new PixiCache(config)

        // Create a test TOML file
        def tomlFile = tempDir.resolve('test-env.toml')
        tomlFile.text = '''
[project]
name = "test-integration"
version = "0.1.0"
description = "Integration test environment"
channels = ["conda-forge"]
platforms = ["linux-64", "osx-64", "osx-arm64"]

[dependencies]
python = ">=3.9"
'''.stripIndent()

        when:
        def prefix = cache.pixiPrefixPath(tomlFile.toString())

        then:
        prefix != null
        prefix.parent == tempDir
        prefix.fileName.toString().startsWith('test-env-')

        cleanup:
        tempDir?.deleteDir()
    }

    @IgnoreIf({ !hasPixiInstalled() })
    def 'should handle pixi environment creation with custom cache directory'() {
        given:
        def customCacheDir = Files.createTempDirectory('custom-pixi-cache')
        def config = new PixiConfig([
            cacheDir: customCacheDir.toString(),
            createTimeout: '10min',
            createOptions: '--no-lockfile-update'
        ], [:])
        def cache = new PixiCache(config)

        when:
        def cacheDir = cache.getCacheDir()

        then:
        cacheDir == customCacheDir
        cacheDir.exists()
        cache.createTimeout == Duration.of('10min')
        cache.createOptions == '--no-lockfile-update'

        cleanup:
        customCacheDir?.deleteDir()
    }

    @IgnoreIf({ !hasPixiInstalled() })
    def 'should validate TOML file detection'() {
        given:
        def cache = new PixiCache(new PixiConfig([:], [:]))

        expect:
        cache.isTomlFilePath('pixi.toml')
        cache.isTomlFilePath('pyproject.toml')
        cache.isTomlFilePath('/path/to/environment.toml')
        !cache.isTomlFilePath('python>=3.9')
        !cache.isTomlFilePath('environment.yaml')
        !cache.isTomlFilePath('multiline\nstring')
    }

    @IgnoreIf({ !hasPixiInstalled() })
    def 'should validate lock file detection'() {
        given:
        def cache = new PixiCache(new PixiConfig([:], [:]))

        expect:
        cache.isLockFilePath('pixi.lock')
        cache.isLockFilePath('/path/to/environment.lock')
        !cache.isLockFilePath('python>=3.9')
        !cache.isLockFilePath('environment.toml')
        !cache.isLockFilePath('multiline\nstring')
    }

    @IgnoreIf({ !hasPixiInstalled() })
    def 'should handle existing prefix directory'() {
        given:
        def tempDir = Files.createTempDirectory('pixi-prefix-test')
        def config = new PixiConfig([:], [:])
        def cache = new PixiCache(config)

        when:
        def prefix = cache.pixiPrefixPath(tempDir.toString())

        then:
        prefix == tempDir
        prefix.isDirectory()

        cleanup:
        tempDir?.deleteDir()
    }

    @IgnoreIf({ !hasPixiInstalled() })
    def 'should handle environment variable cache directory'() {
        given:
        def envCacheDir = Files.createTempDirectory('env-pixi-cache')
        def config = new PixiConfig([:], [NXF_PIXI_CACHEDIR: envCacheDir.toString()])
        def cache = Spy(PixiCache, constructorArgs: [config]) {
            getEnv() >> [NXF_PIXI_CACHEDIR: envCacheDir.toString()]
        }

        when:
        def cacheDir = cache.getCacheDir()

        then:
        cacheDir == envCacheDir
        cacheDir.exists()

        cleanup:
        envCacheDir?.deleteDir()
    }

    @IgnoreIf({ !hasPixiInstalled() })
    def 'should create cache from lock file'() {
        given:
        def tempDir = Files.createTempDirectory('pixi-lock-cache-test')
        def config = new PixiConfig([cacheDir: tempDir.toString()], [:])
        def cache = new PixiCache(config)

        // Create a test lock file (simplified format)
        def lockFile = tempDir.resolve('test.lock')
        lockFile.text = '''
version: 4
environments:
  default:
    channels:
    - url: https://conda.anaconda.org/conda-forge/
    packages:
      linux-64:
      - conda: https://conda.anaconda.org/conda-forge/linux-64/python-3.11.6-hab00c5b_0_cpython.conda
'''.stripIndent()

        when:
        def prefix = cache.pixiPrefixPath(lockFile.toString())

        then:
        prefix != null
        prefix.parent == tempDir
        prefix.fileName.toString().startsWith('test-')

        cleanup:
        tempDir?.deleteDir()
    }

    @IgnoreIf({ !hasPixiInstalled() })
    def 'should reject invalid environment specifications'() {
        given:
        def config = new PixiConfig([:], [:])
        def cache = new PixiCache(config)

        when:
        cache.pixiPrefixPath('invalid\nmultiline\nspec')

        then:
        thrown(IllegalArgumentException)
    }

    @IgnoreIf({ !hasPixiInstalled() })
    def 'should reject non-existent TOML file'() {
        given:
        def config = new PixiConfig([:], [:])
        def cache = new PixiCache(config)

        when:
        cache.pixiPrefixPath('/non/existent/file.toml')

        then:
        thrown(IllegalArgumentException)
    }

    @IgnoreIf({ !hasPixiInstalled() })
    def 'should reject non-existent lock file'() {
        given:
        def config = new PixiConfig([:], [:])
        def cache = new PixiCache(config)

        when:
        cache.pixiPrefixPath('/non/existent/file.lock')

        then:
        thrown(IllegalArgumentException)
    }

    @IgnoreIf({ !hasPixiInstalled() })
    def 'should handle complex TOML file with multiple dependencies'() {
        given:
        def tempDir = Files.createTempDirectory('pixi-complex-toml-test')
        def config = new PixiConfig([cacheDir: tempDir.toString()], [:])
        def cache = new PixiCache(config)

        // Create a complex TOML file
        def tomlFile = tempDir.resolve('complex-env.toml')
        tomlFile.text = '''
[project]
name = "complex-integration-test"
version = "1.0.0"
description = "Complex integration test environment with multiple dependencies"
channels = ["conda-forge", "bioconda", "pytorch"]
platforms = ["linux-64", "osx-64", "osx-arm64"]

[dependencies]
python = ">=3.9,<3.12"
numpy = ">=1.20"
pandas = ">=1.3"
matplotlib = ">=3.5"
scipy = ">=1.7"

[pypi-dependencies]
requests = ">=2.25"

[tasks]
test = "python -c 'import numpy, pandas, matplotlib, scipy; print(\"All packages imported successfully\")'"

[feature.cuda.dependencies]
pytorch = { version = ">=1.12", channel = "pytorch" }
'''.stripIndent()

        when:
        def prefix = cache.pixiPrefixPath(tomlFile.toString())

        then:
        prefix != null
        prefix.parent == tempDir
        prefix.fileName.toString().startsWith('complex-env-')

        cleanup:
        tempDir?.deleteDir()
    }
}
