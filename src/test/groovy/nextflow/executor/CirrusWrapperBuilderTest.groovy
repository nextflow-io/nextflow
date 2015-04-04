/*
 * Copyright (c) 2013-2015, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2015, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
 */

package nextflow.executor

import java.nio.file.Files

import nextflow.processor.TaskConfig
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CirrusWrapperBuilderTest extends Specification {

    def 'test cirrus script wrapper' () {

        given:
        def folder = Files.createTempDirectory('test')
        def processor = Mock(TaskProcessor)
        processor.getProcessEnvironment() >> [:]
        def config = new TaskConfig()
        def task = new TaskRun(workDir: folder, processor: processor, config: config, script: 'echo Hello world!')

        /*
         * simple bash run
         */
        when:
        def bash = new CirrusWrapperBuilder(task)
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
                # fetch scripts
                es3 test s3:/${folder}/.command.env && es3 cat s3:/${folder}/.command.env > .command.env
                es3 test s3:/${folder}/.command.sh && es3 cat s3:/${folder}/.command.sh > .command.sh
                es3 test s3:/${folder}/.command.in && es3 cat s3:/${folder}/.command.in > .command.in
                es3 test s3:/${folder}/.command.run.1 && es3 cat s3:/${folder}/.command.run.1 > .command.run.1

                set -e
                set -u
                NXF_DEBUG=\${NXF_DEBUG:=0}; [[ \$NXF_DEBUG > 2 ]] && set -x

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

                on_exit() {
                  exit_status=\${ret:=\$?}
                  printf \$exit_status > .exitcode && es3 -q -v 0 --no-stats sync .exitcode s3:/${folder} || true
                  sync
                  exit \$exit_status
                }

                on_term() {
                    set +e
                    [[ "\$pid" ]] && nxf_kill \$pid
                }

                trap on_exit EXIT
                trap on_term TERM INT USR1 USR2

                [[ \$NXF_DEBUG > 0 ]] && nxf_env
                es3 touch s3:/${folder}/.command.begin
                [ -f .command.env ] && source .command.env
                NXF_SCRATCH=\${TMPDIR:-`mktemp -d`} && cd \$NXF_SCRATCH

                set +e
                (
                /bin/bash -ue .command.sh
                ) > >(tee .command.out) 2> >(tee .command.err >&2) &
                pid=\$!
                wait \$pid || ret=\$?
                es3 -q -v 0 --no-stats sync .command.out s3:/${folder} || true
                es3 -q -v 0 --no-stats sync .command.err s3:/${folder} || true
                """
                        .stripIndent().leftTrim()


        cleanup:
        folder?.deleteDir()
    }

}
