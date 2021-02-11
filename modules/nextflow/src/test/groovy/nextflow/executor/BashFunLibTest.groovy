package nextflow.executor

import java.util.concurrent.TimeUnit
import java.nio.file.Files

import nextflow.util.Duration
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class BashFunLibTest extends Specification {

    def 'should fail on errors in nxf_parallel' () {

        given:
        def scriptFile = Files.createTempFile("test", "sh")
        def script = """
        #!/bin/bash

        set -e
        """.stripIndent() +

        BashFunLib.body(5, 5, Duration.of('1 sec')) +

        """
        cmds=()
        cmds+=("true")
        cmds+=("false")
        nxf_parallel "\${cmds[@]}"
        """.stripIndent()

        scriptFile.text = script

        def process = "bash ${scriptFile}".execute()
        process.waitFor(1, TimeUnit.SECONDS)

        expect:
        process.exitValue() == 1

        cleanup:
        if( scriptFile ) Files.delete(scriptFile)

    }

    def 'should succeed with nxf_parallel' () {

        given:
        def scriptFile = Files.createTempFile("test", "sh")
        def script = """
        #!/bin/bash

        set -e
        """.stripIndent() +

        BashFunLib.body(5, 5, Duration.of('1 sec')) +

        """
        cmds=()
        cmds+=("true")
        cmds+=("true")
        nxf_parallel "\${cmds[@]}"
        """.stripIndent()

        scriptFile.text = script

        def process = "bash ${scriptFile}".execute()
        process.waitFor(1, TimeUnit.SECONDS)

        expect:
        process.exitValue() == 0

        cleanup:
        if( scriptFile ) Files.delete(scriptFile)

    }

}
