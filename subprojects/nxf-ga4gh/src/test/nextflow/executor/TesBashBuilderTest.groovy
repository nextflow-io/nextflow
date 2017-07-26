package nextflow.executor

import java.nio.file.Files

import nextflow.processor.TaskBean
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TesBashBuilderTest extends Specification {


    def 'test bash wrapper' () {

        given:
        def folder = Files.createTempDirectory('test')

        /*
         * simple bash run
         */
        when:
        def bash = new TesBashBuilder([
                name: 'Hello 1',
                workDir: folder,
                script: 'echo Hello world!',
                headerScript: '#BSUB -x 1\n#BSUB -y 2'
        ] as TaskBean )
        bash.build()

        then:
        Files.exists(folder.resolve('.command.sh'))
        Files.exists(folder.resolve('.command.run'))

        folder.resolve('.command.sh').text ==
                '''
                #!/bin/bash -ue
                echo Hello world!
                '''
                        .stripIndent().leftTrim()


        folder.resolve('.command.run').text ==
                """
                #!/bin/bash
                #BSUB -x 1
                #BSUB -y 2
                # NEXTFLOW TASK: Hello 1
                set -e
                set -u
                NXF_DEBUG=\${NXF_DEBUG:=0}; [[ \$NXF_DEBUG > 1 ]] && set -x

                nxf_env() {
                    echo '============= task environment ============='
                    env | sort | sed "s/\\(.*\\)AWS\\(.*\\)=\\(.\\{6\\}\\).*/\\1AWS\\2=\\3xxxxxxxxxxxxx/"
                    echo '============= task output =================='
                }

                nxf_kill() {
                    declare -a ALL_CHILD
                    while read P PP;do
                        ALL_CHILD[\$PP]+=" \$P"
                    done < <(ps -e -o pid= -o ppid=)

                    walk() {
                        [[ \$1 != \$\$ ]] && kill \$1 2>/dev/null || true
                        for i in \${ALL_CHILD[\$1]:=}; do walk \$i; done
                    }

                    walk \$1
                }

                nxf_mktemp() {
                    local base=\${1:-/tmp}
                    if [[ \$(uname) = Darwin ]]; then mktemp -d \$base/nxf.XXXXXXXXXX
                    else TMPDIR="\$base" mktemp -d -t nxf.XXXXXXXXXX
                    fi
                }

                on_exit() {
                  exit_status=\${ret:=\$?}
                  printf \$exit_status > .exitcode
                  set +u
                  [[ "\$COUT" ]] && rm -f "\$COUT" || true
                  [[ "\$CERR" ]] && rm -f "\$CERR" || true
                  exit \$exit_status
                }

                on_term() {
                    set +e
                    [[ "\$pid" ]] && nxf_kill \$pid
                }

                trap on_exit EXIT
                trap on_term TERM INT USR1 USR2

                NXF_SCRATCH=''
                [[ \$NXF_DEBUG > 0 ]] && nxf_env

                [ -f .command.env ] && source .command.env
                [[ \$NXF_SCRATCH ]] && echo "nxf-scratch-dir \$HOSTNAME:\$NXF_SCRATCH" && cd \$NXF_SCRATCH

                set +e
                COUT=\$PWD/.command.po; mkfifo "\$COUT"
                CERR=\$PWD/.command.pe; mkfifo "\$CERR"
                tee .command.out < "\$COUT" &
                tee1=\$!
                tee .command.err < "\$CERR" >&2 &
                tee2=\$!
                (
                /bin/bash -ue .command.sh
                ) >"\$COUT" 2>"\$CERR" &
                pid=\$!
                wait \$pid || ret=\$?
                wait \$tee1 \$tee2
                """
                        .stripIndent().leftTrim()


        cleanup:
        folder?.deleteDir()
    }


}
