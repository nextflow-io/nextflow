package nextflow.executor

import nextflow.util.Duration

/**
 * Bash common functions library
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class BashFunLib {

    static String body(int maxConnect, int maxAttempts, Duration delayBetweenAttempts) {
        """\
        # bash helper functions
        nxf_cp_retry() {
            local max_attempts=$maxAttempts
            local timeout=${delayBetweenAttempts.seconds}
            local attempt=0
            local exitCode=0
            while (( \$attempt < \$max_attempts ))
            do
              if "\$@"
                then
                  return 0
              else
                exitCode=\$?
              fi
              if [[ \$exitCode == 0 ]]
              then
                break
              fi
              sleep \$timeout
              attempt=\$(( attempt + 1 ))
              timeout=\$(( timeout * 2 ))
            done
        }
        
        nxf_parallel() {
            IFS=\$'\\n'
            local cmd=("\$@")
            local cpus=\$(nproc 2>/dev/null || < /proc/cpuinfo grep '^process' -c)
            local max=\$(if (( cpus>$maxConnect )); then echo $maxConnect; else echo \$cpus; fi)
            local i=0
            local pid=()
            (
            set +u
            while ((i<\${#cmd[@]})); do
                local copy=()
                for x in "\${pid[@]}"; do
                  [[ -e /proc/\$x ]] && copy+=(\$x)
                done
                pid=("\${copy[@]}")
                
                if ((\${#pid[@]}>=\$max)); then
                  sleep 1
                else
                  eval "\${cmd[\$i]}" &
                  pid+=(\$!)
                  ((i+=1))
                fi
            done
            ((\${#pid[@]}>0)) && wait \${pid[@]}
            )
            unset IFS
        }
        """.stripIndent()
    }
}
