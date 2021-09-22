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

    def 'should render bash body' () {

        expect:
        BashFunLib.body(123, 321, Duration.of('1 sec')) == '''\
            # bash helper functions
            nxf_cp_retry() {
                local max_attempts=321
                local timeout=1
                local attempt=0
                local exitCode=0
                while (( $attempt < $max_attempts ))
                do
                  if "$@"
                    then
                      return 0
                  else
                    exitCode=$?
                  fi
                  if [[ $exitCode == 0 ]]
                  then
                    break
                  fi
                  nxf_sleep $timeout
                  attempt=$(( attempt + 1 ))
                  timeout=$(( timeout * 2 ))
                done
            }
            
            nxf_parallel() {
                IFS=$'\\n'
                local cmd=("$@")
                local cpus=$(nproc 2>/dev/null || < /proc/cpuinfo grep '^process' -c)
                local max=$(if (( cpus>123 )); then echo 123; else echo $cpus; fi)
                local i=0
                local pid=()
                (
                set +u
                while ((i<${#cmd[@]})); do
                    local copy=()
                    for x in "${pid[@]}"; do
                      [[ -e /proc/$x ]] && copy+=($x)
                    done
                    pid=("${copy[@]}")
            
                    if ((${#pid[@]}>=$max)); then
                      nxf_sleep 0.2
                    else
                      eval "${cmd[$i]}" &
                      pid+=($!)
                      ((i+=1))
                    fi
                done
                for p in "${pid[@]}"; do
                    wait $p
                done
                )
                unset IFS
            }
            '''.stripIndent()
    }


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
        process.waitFor(30, TimeUnit.SECONDS)

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
        process.waitFor(30, TimeUnit.SECONDS)

        expect:
        process.exitValue() == 0

        cleanup:
        if( scriptFile ) Files.delete(scriptFile)

    }

}
