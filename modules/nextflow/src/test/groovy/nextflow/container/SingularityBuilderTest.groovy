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

package nextflow.container

import java.nio.file.Paths

import nextflow.SysEnv
import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class SingularityBuilderTest extends Specification {

    def 'should get the run command line and auto mounts disable' () {
        given:
        SysEnv.push(NXF_SINGULARITY_RUN_COMMAND:'run', NXF_SINGULARITY_AUTO_MOUNTS:'false')
        and:
        def path1 = Paths.get('/foo/data/file1')
        def path2 = Paths.get('/bar/data/file2')
        def path3 = Paths.get('/bar/data file')

        expect:
        new SingularityBuilder('busybox')
                .build()
                .runCommand == 'set +u; env - PATH="$PATH" ${TMP:+SINGULARITYENV_TMP="$TMP"} ${TMPDIR:+SINGULARITYENV_TMPDIR="$TMPDIR"} singularity run --no-home --pid busybox'

        new SingularityBuilder('busybox')
                .params(engineOptions: '-q -v')
                .build()
                .runCommand == 'set +u; env - PATH="$PATH" ${TMP:+SINGULARITYENV_TMP="$TMP"} ${TMPDIR:+SINGULARITYENV_TMPDIR="$TMPDIR"} singularity -q -v run --no-home --pid busybox'

        new SingularityBuilder('busybox')
                .params(runOptions: '--contain --writable')
                .build()
                .runCommand == 'set +u; env - PATH="$PATH" ${TMP:+SINGULARITYENV_TMP="$TMP"} ${TMPDIR:+SINGULARITYENV_TMPDIR="$TMPDIR"} singularity run --no-home --pid --contain --writable busybox'

        new SingularityBuilder('ubuntu')
                .addMount(path1)
                .build()
                .runCommand == 'set +u; env - PATH="$PATH" ${TMP:+SINGULARITYENV_TMP="$TMP"} ${TMPDIR:+SINGULARITYENV_TMPDIR="$TMPDIR"} singularity run --no-home --pid ubuntu'

        new SingularityBuilder('ubuntu')
                .addMount(path1)
                .addMount(path2)
                .params(autoMounts: true)
                .build()
                .runCommand == 'set +u; env - PATH="$PATH" ${TMP:+SINGULARITYENV_TMP="$TMP"} ${TMPDIR:+SINGULARITYENV_TMPDIR="$TMPDIR"} singularity run --no-home --pid -B /foo/data/file1 -B /bar/data/file2 -B "$NXF_TASK_WORKDIR" ubuntu'

        new SingularityBuilder('ubuntu')
                .addMount(path1)
                .addMount(path1)
                .params(autoMounts: true)
                .build()
                .runCommand == 'set +u; env - PATH="$PATH" ${TMP:+SINGULARITYENV_TMP="$TMP"} ${TMPDIR:+SINGULARITYENV_TMPDIR="$TMPDIR"} singularity run --no-home --pid -B /foo/data/file1 -B "$NXF_TASK_WORKDIR" ubuntu'

        new SingularityBuilder('ubuntu')
                .addMount(path1)
                .addMount(path1)
                .params(autoMounts: true)
                .params(readOnlyInputs: true)
                .build()
                .runCommand == 'set +u; env - PATH="$PATH" ${TMP:+SINGULARITYENV_TMP="$TMP"} ${TMPDIR:+SINGULARITYENV_TMPDIR="$TMPDIR"} singularity run --no-home --pid -B /foo/data/file1:/foo/data/file1:ro -B "$NXF_TASK_WORKDIR" ubuntu'

        new SingularityBuilder('ubuntu')
                .addMount(path3)
                .params(autoMounts: true)
                .build()
                .runCommand == 'set +u; env - PATH="$PATH" ${TMP:+SINGULARITYENV_TMP="$TMP"} ${TMPDIR:+SINGULARITYENV_TMPDIR="$TMPDIR"} singularity run --no-home --pid -B /bar/data\\ file -B "$NXF_TASK_WORKDIR" ubuntu'

        new SingularityBuilder('ubuntu')
                .params(newPidNamespace: false)
                .build()
                .runCommand == 'set +u; env - PATH="$PATH" ${TMP:+SINGULARITYENV_TMP="$TMP"} ${TMPDIR:+SINGULARITYENV_TMPDIR="$TMPDIR"} singularity run --no-home ubuntu'

        cleanup:
        SysEnv.pop()
    }

    def 'should get the run command line' () {

        given:
        def path1 = Paths.get('/foo/data/file1')
        def path2 = Paths.get('/bar/data/file2')
        def path3 = Paths.get('/bar/data file')

        expect:
        new SingularityBuilder('busybox')
                .build()
                .runCommand == 'set +u; env - PATH="$PATH" ${TMP:+SINGULARITYENV_TMP="$TMP"} ${TMPDIR:+SINGULARITYENV_TMPDIR="$TMPDIR"} singularity exec --no-home --pid -B "$NXF_TASK_WORKDIR" busybox'

        new SingularityBuilder('busybox')
                .params(engineOptions: '-q -v')
                .build()
                .runCommand == 'set +u; env - PATH="$PATH" ${TMP:+SINGULARITYENV_TMP="$TMP"} ${TMPDIR:+SINGULARITYENV_TMPDIR="$TMPDIR"} singularity -q -v exec --no-home --pid -B "$NXF_TASK_WORKDIR" busybox'

        new SingularityBuilder('busybox')
                .params(runOptions: '--contain --writable')
                .build()
                .runCommand == 'set +u; env - PATH="$PATH" ${TMP:+SINGULARITYENV_TMP="$TMP"} ${TMPDIR:+SINGULARITYENV_TMPDIR="$TMPDIR"} singularity exec --no-home --pid -B "$NXF_TASK_WORKDIR" --contain --writable busybox'

        new SingularityBuilder('ubuntu')
                .params(autoMounts: false)
                .build()
                .runCommand == 'set +u; env - PATH="$PATH" ${TMP:+SINGULARITYENV_TMP="$TMP"} ${TMPDIR:+SINGULARITYENV_TMPDIR="$TMPDIR"} singularity exec --no-home --pid ubuntu'

        new SingularityBuilder('ubuntu')
                .addMount(path1)
                .addMount(path2)
                .params(autoMounts: true)
                .build()
                .runCommand == 'set +u; env - PATH="$PATH" ${TMP:+SINGULARITYENV_TMP="$TMP"} ${TMPDIR:+SINGULARITYENV_TMPDIR="$TMPDIR"} singularity exec --no-home --pid -B /foo/data/file1 -B /bar/data/file2 -B "$NXF_TASK_WORKDIR" ubuntu'

        new SingularityBuilder('ubuntu')
                .addMount(path1)
                .addMount(path1)
                .params(autoMounts: true)
                .build()
                .runCommand == 'set +u; env - PATH="$PATH" ${TMP:+SINGULARITYENV_TMP="$TMP"} ${TMPDIR:+SINGULARITYENV_TMPDIR="$TMPDIR"} singularity exec --no-home --pid -B /foo/data/file1 -B "$NXF_TASK_WORKDIR" ubuntu'

        new SingularityBuilder('ubuntu')
                .addMount(path1)
                .addMount(path1)
                .params(autoMounts: true)
                .params(readOnlyInputs: true)
                .build()
                .runCommand == 'set +u; env - PATH="$PATH" ${TMP:+SINGULARITYENV_TMP="$TMP"} ${TMPDIR:+SINGULARITYENV_TMPDIR="$TMPDIR"} singularity exec --no-home --pid -B /foo/data/file1:/foo/data/file1:ro -B "$NXF_TASK_WORKDIR" ubuntu'

        new SingularityBuilder('ubuntu')
                .addMount(path3)
                .params(autoMounts: true)
                .build()
                .runCommand == 'set +u; env - PATH="$PATH" ${TMP:+SINGULARITYENV_TMP="$TMP"} ${TMPDIR:+SINGULARITYENV_TMPDIR="$TMPDIR"} singularity exec --no-home --pid -B /bar/data\\ file -B "$NXF_TASK_WORKDIR" ubuntu'

        new SingularityBuilder('ubuntu')
                .params(newPidNamespace: false)
                .params(autoMounts: false)
                .build()
                .runCommand == 'set +u; env - PATH="$PATH" ${TMP:+SINGULARITYENV_TMP="$TMP"} ${TMPDIR:+SINGULARITYENV_TMPDIR="$TMPDIR"} singularity exec --no-home ubuntu'

        new SingularityBuilder('ubuntu')
            .params(oci: true)
            .build()
            .runCommand == 'set +u; env - PATH="$PATH" ${TMP:+SINGULARITYENV_TMP="$TMP"} ${TMPDIR:+SINGULARITYENV_TMPDIR="$TMPDIR"} ${XDG_RUNTIME_DIR:+XDG_RUNTIME_DIR="$XDG_RUNTIME_DIR"} ${DBUS_SESSION_BUS_ADDRESS:+DBUS_SESSION_BUS_ADDRESS="$DBUS_SESSION_BUS_ADDRESS"} singularity exec --no-home --oci -B "$NXF_TASK_WORKDIR" ubuntu'

        new SingularityBuilder('ubuntu')
            .params(ociMode: true)
            .build()
            .runCommand == 'set +u; env - PATH="$PATH" ${TMP:+SINGULARITYENV_TMP="$TMP"} ${TMPDIR:+SINGULARITYENV_TMPDIR="$TMPDIR"} ${XDG_RUNTIME_DIR:+XDG_RUNTIME_DIR="$XDG_RUNTIME_DIR"} ${DBUS_SESSION_BUS_ADDRESS:+DBUS_SESSION_BUS_ADDRESS="$DBUS_SESSION_BUS_ADDRESS"} singularity exec --no-home --oci -B "$NXF_TASK_WORKDIR" ubuntu'


    }

    def 'should mount home directory if specified' () {
        when:
        SysEnv.push(NXF_SINGULARITY_HOME_MOUNT: FLAG)

        then:
        new SingularityBuilder('ubuntu')
                .build()
                .runCommand
                .contains('--no-home ') == OPTS

        cleanup:
        SysEnv.pop()

        where:
        FLAG << ['false', 'true']
        OPTS << [true, false]
    }

    def 'should set new pid namespace' () {
        when:
        SysEnv.push(NXF_SINGULARITY_NEW_PID_NAMESPACE: FLAG)

        then:
        new SingularityBuilder('ubuntu')
                .build()
                .runCommand
                .contains('--pid ') == OPTS

        cleanup:
        SysEnv.pop()

        where:
        FLAG << ['false', 'true']
        OPTS << [false, true]
    }

    def 'should return export variables' () {

        expect:
        new SingularityBuilder('busybox')
                .getEnvExports() == ''

        new SingularityBuilder('busybox')
                .addEnv('X=1')
                .addEnv(ALPHA:'aaa', BETA: 'bbb')
                .build()
                .runCommand == 'set +u; env - PATH="$PATH" ${TMP:+SINGULARITYENV_TMP="$TMP"} ${TMPDIR:+SINGULARITYENV_TMPDIR="$TMPDIR"} SINGULARITYENV_X="1" SINGULARITYENV_ALPHA="aaa" SINGULARITYENV_BETA="bbb" singularity exec --no-home --pid -B "$NXF_TASK_WORKDIR" busybox'

        new SingularityBuilder('busybox')
                .addEnv('CUDA_VISIBLE_DEVICES')
                .build()
                .runCommand == 'set +u; env - PATH="$PATH" ${TMP:+SINGULARITYENV_TMP="$TMP"} ${TMPDIR:+SINGULARITYENV_TMPDIR="$TMPDIR"} ${CUDA_VISIBLE_DEVICES:+SINGULARITYENV_CUDA_VISIBLE_DEVICES="$CUDA_VISIBLE_DEVICES"} singularity exec --no-home --pid -B "$NXF_TASK_WORKDIR" busybox'

    }


    def 'should get run command' () {

        when:
        def cmd = new SingularityBuilder('ubuntu.img').build().getRunCommand()
        then:
        cmd == 'set +u; env - PATH="$PATH" ${TMP:+SINGULARITYENV_TMP="$TMP"} ${TMPDIR:+SINGULARITYENV_TMPDIR="$TMPDIR"} singularity exec --no-home --pid -B "$NXF_TASK_WORKDIR" ubuntu.img'

        when:
        cmd = new SingularityBuilder('ubuntu.img').build().getRunCommand('bwa --this --that file.fastq')
        then:
        cmd == 'set +u; env - PATH="$PATH" ${TMP:+SINGULARITYENV_TMP="$TMP"} ${TMPDIR:+SINGULARITYENV_TMPDIR="$TMPDIR"} singularity exec --no-home --pid -B "$NXF_TASK_WORKDIR" ubuntu.img bwa --this --that file.fastq'

        when:
        cmd = new SingularityBuilder('ubuntu.img').params(entry:'/bin/sh').build().getRunCommand('bwa --this --that file.fastq')
        then:
        cmd == 'set +u; env - PATH="$PATH" ${TMP:+SINGULARITYENV_TMP="$TMP"} ${TMPDIR:+SINGULARITYENV_TMPDIR="$TMPDIR"} singularity exec --no-home --pid -B "$NXF_TASK_WORKDIR" ubuntu.img /bin/sh -c "cd $NXF_TASK_WORKDIR; bwa --this --that file.fastq"'
    }

    @Unroll
    def 'test singularity env'() {

        given:
        def builder = Spy(SingularityBuilder)

        expect:
        builder.makeEnv(ENV).toString() == RESULT

        where:
        ENV                         | RESULT
        'X=1'                       | 'SINGULARITYENV_X="1"'
        'BAR'                       | '${BAR:+SINGULARITYENV_BAR="$BAR"}'
        [VAR_X:1, VAR_Y: 2]         | 'SINGULARITYENV_VAR_X="1" SINGULARITYENV_VAR_Y="2"'
        [SINGULARITY_BIND: 'foo', SINGULARITYENV_FOO: 'x', BAR: 'y'] | 'SINGULARITY_BIND="foo" SINGULARITYENV_FOO="x" SINGULARITYENV_BAR="y"'
        'SINGULARITY_FOO'          | '${SINGULARITY_FOO:+SINGULARITY_FOO="$SINGULARITY_FOO"}'
        'SINGULARITYENV_FOO'       | '${SINGULARITYENV_FOO:+SINGULARITYENV_FOO="$SINGULARITYENV_FOO"}'
        'SINGULARITYENV_X=1'       | 'SINGULARITYENV_X="1"'
    }

    @Unroll
    def 'should quote env value' () {
        given:
        def builder = Spy(SingularityBuilder)

        expect:
        builder.quoteValue(STR) == EXPECTED

        where:
        STR             | EXPECTED
        null            | null
        ''              | ''
        'foo'           | '"foo"'
        '"foo"'         | '"foo"'
        "'foo'"         | "'foo'"
        and:
        'X=foo'         | 'X="foo"'
        'X="foo"'       | 'X="foo"'
        'X='            | 'X='
    }
}
