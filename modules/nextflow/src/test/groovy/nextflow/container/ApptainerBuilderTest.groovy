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

package nextflow.container

import java.nio.file.Paths

import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ApptainerBuilderTest extends Specification {

    def 'should get the exec command line' () {

        given:
        def path1 = Paths.get('/foo/data/file1')
        def path2 = Paths.get('/bar/data/file2')
        def path3 = Paths.get('/bar/data file')

        expect:
        new ApptainerBuilder('busybox')
                .build()
                .runCommand == 'set +u; env - PATH="$PATH" ${TMP:+APPTAINERENV_TMP="$TMP"} ${TMPDIR:+APPTAINERENV_TMPDIR="$TMPDIR"} apptainer exec busybox'

        new ApptainerBuilder('busybox')
                .params(engineOptions: '-q -v')
                .build()
                .runCommand == 'set +u; env - PATH="$PATH" ${TMP:+APPTAINERENV_TMP="$TMP"} ${TMPDIR:+APPTAINERENV_TMPDIR="$TMPDIR"} apptainer -q -v exec busybox'

        new ApptainerBuilder('busybox')
                .params(runOptions: '--contain --writable')
                .build()
                .runCommand == 'set +u; env - PATH="$PATH" ${TMP:+APPTAINERENV_TMP="$TMP"} ${TMPDIR:+APPTAINERENV_TMPDIR="$TMPDIR"} apptainer exec --contain --writable busybox'

        new ApptainerBuilder('ubuntu')
                .addMount(path1)
                .build()
                .runCommand == 'set +u; env - PATH="$PATH" ${TMP:+APPTAINERENV_TMP="$TMP"} ${TMPDIR:+APPTAINERENV_TMPDIR="$TMPDIR"} apptainer exec ubuntu'

        new ApptainerBuilder('ubuntu')
                .addMount(path1)
                .addMount(path2)
                .params(autoMounts: true)
                .build()
                .runCommand == 'set +u; env - PATH="$PATH" ${TMP:+APPTAINERENV_TMP="$TMP"} ${TMPDIR:+APPTAINERENV_TMPDIR="$TMPDIR"} apptainer exec -B /foo/data/file1 -B /bar/data/file2 -B "$PWD" ubuntu'

        new ApptainerBuilder('ubuntu')
                .addMount(path1)
                .addMount(path1)
                .params(autoMounts: true)
                .build()
                .runCommand == 'set +u; env - PATH="$PATH" ${TMP:+APPTAINERENV_TMP="$TMP"} ${TMPDIR:+APPTAINERENV_TMPDIR="$TMPDIR"} apptainer exec -B /foo/data/file1 -B "$PWD" ubuntu'

        new ApptainerBuilder('ubuntu')
                .addMount(path1)
                .addMount(path1)
                .params(autoMounts: true)
                .params(readOnlyInputs: true)
                .build()
                .runCommand == 'set +u; env - PATH="$PATH" ${TMP:+APPTAINERENV_TMP="$TMP"} ${TMPDIR:+APPTAINERENV_TMPDIR="$TMPDIR"} apptainer exec -B /foo/data/file1:/foo/data/file1:ro -B "$PWD" ubuntu'

        new ApptainerBuilder('ubuntu')
                .addMount(path3)
                .params(autoMounts: true)
                .build()
                .runCommand == 'set +u; env - PATH="$PATH" ${TMP:+APPTAINERENV_TMP="$TMP"} ${TMPDIR:+APPTAINERENV_TMPDIR="$TMPDIR"} apptainer exec -B /bar/data\\ file -B "$PWD" ubuntu'

    }

    def 'should return export variables' () {

        expect:
        new ApptainerBuilder('busybox')
                .getEnvExports() == ''

        new ApptainerBuilder('busybox')
                .addEnv('X=1')
                .addEnv(ALPHA:'aaa', BETA: 'bbb')
                .build()
                .runCommand == 'set +u; env - PATH="$PATH" ${TMP:+APPTAINERENV_TMP="$TMP"} ${TMPDIR:+APPTAINERENV_TMPDIR="$TMPDIR"} APPTAINERENV_X=1 APPTAINERENV_ALPHA="aaa" APPTAINERENV_BETA="bbb" apptainer exec busybox'

        new ApptainerBuilder('busybox')
                .addEnv('CUDA_VISIBLE_DEVICES')
                .build()
                .runCommand == 'set +u; env - PATH="$PATH" ${TMP:+APPTAINERENV_TMP="$TMP"} ${TMPDIR:+APPTAINERENV_TMPDIR="$TMPDIR"} ${CUDA_VISIBLE_DEVICES:+APPTAINERENV_CUDA_VISIBLE_DEVICES="$CUDA_VISIBLE_DEVICES"} apptainer exec busybox'

    }


    def 'should get run command' () {

        when:
        def cmd = new ApptainerBuilder('ubuntu.img').build().getRunCommand()
        then:
        cmd == 'set +u; env - PATH="$PATH" ${TMP:+APPTAINERENV_TMP="$TMP"} ${TMPDIR:+APPTAINERENV_TMPDIR="$TMPDIR"} apptainer exec ubuntu.img'

        when:
        cmd = new ApptainerBuilder('ubuntu.img').build().getRunCommand('bwa --this --that file.fastq')
        then:
        cmd == 'set +u; env - PATH="$PATH" ${TMP:+APPTAINERENV_TMP="$TMP"} ${TMPDIR:+APPTAINERENV_TMPDIR="$TMPDIR"} apptainer exec ubuntu.img bwa --this --that file.fastq'

        when:
        cmd = new ApptainerBuilder('ubuntu.img').params(entry:'/bin/sh').build().getRunCommand('bwa --this --that file.fastq')
        then:
        cmd == 'set +u; env - PATH="$PATH" ${TMP:+APPTAINERENV_TMP="$TMP"} ${TMPDIR:+APPTAINERENV_TMPDIR="$TMPDIR"} apptainer exec ubuntu.img /bin/sh -c "cd $PWD; bwa --this --that file.fastq"'
    }

    @Unroll
    def 'test apptainer env'() {

        given:
        def builder = Spy(ApptainerBuilder)

        expect:
        builder.makeEnv(ENV).toString() == RESULT

        where:
        ENV                         | RESULT
        'X=1'                       | 'APPTAINERENV_X=1'
        'BAR'                       | '${BAR:+APPTAINERENV_BAR="$BAR"}'
        [VAR_X:1, VAR_Y: 2]         | 'APPTAINERENV_VAR_X="1" APPTAINERENV_VAR_Y="2"'
        [APPTAINER_BIND: 'foo', APPTAINERENV_FOO: 'x', BAR: 'y'] | 'APPTAINER_BIND="foo" APPTAINERENV_FOO="x" APPTAINERENV_BAR="y"'
        'APPTAINER_FOO'          | '${APPTAINER_FOO:+APPTAINER_FOO="$APPTAINER_FOO"}'
        'APPTAINERENV_FOO'       | '${APPTAINERENV_FOO:+APPTAINERENV_FOO="$APPTAINERENV_FOO"}'
        'APPTAINERENV_X=1'       | 'APPTAINERENV_X=1'
    }
}
