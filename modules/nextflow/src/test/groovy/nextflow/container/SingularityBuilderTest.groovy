/*
 * Copyright 2020, Seqera Labs
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
class SingularityBuilderTest extends Specification {

    def 'should get the exec command line' () {

        given:
        def path1 = Paths.get('/foo/data/file1')
        def path2 = Paths.get('/bar/data/file2')
        def path3 = Paths.get('/bar/data file')

        expect:
        new SingularityBuilder('busybox')
                .build()
                .runCommand == 'set +u; env - PATH="$PATH" SINGULARITYENV_TMP="$TMP" SINGULARITYENV_TMPDIR="$TMPDIR" singularity exec busybox'

        new SingularityBuilder('busybox')
                .params(engineOptions: '-q -v')
                .build()
                .runCommand == 'set +u; env - PATH="$PATH" SINGULARITYENV_TMP="$TMP" SINGULARITYENV_TMPDIR="$TMPDIR" singularity -q -v exec busybox'

        new SingularityBuilder('busybox')
                .params(runOptions: '--contain --writable')
                .build()
                .runCommand == 'set +u; env - PATH="$PATH" SINGULARITYENV_TMP="$TMP" SINGULARITYENV_TMPDIR="$TMPDIR" singularity exec --contain --writable busybox'

        new SingularityBuilder('ubuntu')
                .addMount(path1)
                .build()
                .runCommand == 'set +u; env - PATH="$PATH" SINGULARITYENV_TMP="$TMP" SINGULARITYENV_TMPDIR="$TMPDIR" singularity exec ubuntu'

        new SingularityBuilder('ubuntu')
                .addMount(path1)
                .addMount(path2)
                .params(autoMounts: true)
                .build()
                .runCommand == 'set +u; env - PATH="$PATH" SINGULARITYENV_TMP="$TMP" SINGULARITYENV_TMPDIR="$TMPDIR" singularity exec -B /foo/data/file1 -B /bar/data/file2 -B "$PWD" ubuntu'

        new SingularityBuilder('ubuntu')
                .addMount(path1)
                .addMount(path1)
                .params(autoMounts: true)
                .build()
                .runCommand == 'set +u; env - PATH="$PATH" SINGULARITYENV_TMP="$TMP" SINGULARITYENV_TMPDIR="$TMPDIR" singularity exec -B /foo/data/file1 -B "$PWD" ubuntu'

        new SingularityBuilder('ubuntu')
                .addMount(path1)
                .addMount(path1)
                .params(autoMounts: true)
                .params(readOnlyInputs: true)
                .build()
                .runCommand == 'set +u; env - PATH="$PATH" SINGULARITYENV_TMP="$TMP" SINGULARITYENV_TMPDIR="$TMPDIR" singularity exec -B /foo/data/file1:/foo/data/file1:ro -B "$PWD" ubuntu'

        new SingularityBuilder('ubuntu')
                .addMount(path3)
                .params(autoMounts: true)
                .build()
                .runCommand == 'set +u; env - PATH="$PATH" SINGULARITYENV_TMP="$TMP" SINGULARITYENV_TMPDIR="$TMPDIR" singularity exec -B /bar/data\\ file -B "$PWD" ubuntu'

    }

    def 'should return export variables' () {

        expect:
        new SingularityBuilder('busybox')
                .getEnvExports() == ''

        new SingularityBuilder('busybox')
                .addEnv('X=1')
                .addEnv(ALPHA:'aaa', BETA: 'bbb')
                .build()
                .runCommand == 'set +u; env - PATH="$PATH" SINGULARITYENV_TMP="$TMP" SINGULARITYENV_TMPDIR="$TMPDIR" SINGULARITYENV_X=1 SINGULARITYENV_ALPHA="aaa" SINGULARITYENV_BETA="bbb" singularity exec busybox'

        new SingularityBuilder('busybox')
                .addEnv('CUDA_VISIBLE_DEVICES')
                .build()
                .runCommand == 'set +u; env - PATH="$PATH" SINGULARITYENV_TMP="$TMP" SINGULARITYENV_TMPDIR="$TMPDIR" ${CUDA_VISIBLE_DEVICES:+SINGULARITYENV_CUDA_VISIBLE_DEVICES="$CUDA_VISIBLE_DEVICES"} singularity exec busybox'

    }


    def 'should get run command' () {

        when:
        def cmd = new SingularityBuilder('ubuntu.img').build().getRunCommand()
        then:
        cmd == 'set +u; env - PATH="$PATH" SINGULARITYENV_TMP="$TMP" SINGULARITYENV_TMPDIR="$TMPDIR" singularity exec ubuntu.img'

        when:
        cmd = new SingularityBuilder('ubuntu.img').build().getRunCommand('bwa --this --that file.fastq')
        then:
        cmd == 'set +u; env - PATH="$PATH" SINGULARITYENV_TMP="$TMP" SINGULARITYENV_TMPDIR="$TMPDIR" singularity exec ubuntu.img bwa --this --that file.fastq'

        when:
        cmd = new SingularityBuilder('ubuntu.img').params(entry:'/bin/sh').build().getRunCommand('bwa --this --that file.fastq')
        then:
        cmd == 'set +u; env - PATH="$PATH" SINGULARITYENV_TMP="$TMP" SINGULARITYENV_TMPDIR="$TMPDIR" singularity exec ubuntu.img /bin/sh -c "cd $PWD; bwa --this --that file.fastq"'
    }


    def 'test memory and cpus'() {

        expect:
        new SingularityBuilder('busybox')
                .setCpus(2)
                .build()
                .runCommand == 'set +u; env - PATH="$PATH" SINGULARITYENV_TMP="$TMP" SINGULARITYENV_TMPDIR="$TMPDIR" singularity exec --vm-cpu 2.0 busybox'


        new SingularityBuilder('busybox')
                .setMemory('100')
                .build()
                .runCommand == 'set +u; env - PATH="$PATH" SINGULARITYENV_TMP="$TMP" SINGULARITYENV_TMPDIR="$TMPDIR" singularity exec --vm-ram 100 busybox'


        new SingularityBuilder('busybox')
                .setCpus(2.5)
                .setMemory('100')
                .build()
                .runCommand == 'set +u; env - PATH="$PATH" SINGULARITYENV_TMP="$TMP" SINGULARITYENV_TMPDIR="$TMPDIR" singularity exec --vm-cpu 2.5 --vm-ram 100 busybox'


    }

    @Unroll
    def 'test singularity env'() {

        given:
        def builder = Spy(SingularityBuilder)

        expect:
        builder.makeEnv(ENV).toString() == RESULT

        where:
        ENV                         | RESULT
        'X=1'                       | 'SINGULARITYENV_X=1'
        'BAR'                       | '${BAR:+SINGULARITYENV_BAR="$BAR"}'
        [VAR_X:1, VAR_Y: 2]         | 'SINGULARITYENV_VAR_X="1" SINGULARITYENV_VAR_Y="2"'
        [SINGULARITY_BIND: 'foo', SINGULARITYENV_FOO: 'x', BAR: 'y'] | 'SINGULARITY_BIND="foo" SINGULARITYENV_FOO="x" SINGULARITYENV_BAR="y"'
        'SINGULARITY_FOO'          | '${SINGULARITY_FOO:+SINGULARITY_FOO="$SINGULARITY_FOO"}'
        'SINGULARITYENV_FOO'       | '${SINGULARITYENV_FOO:+SINGULARITYENV_FOO="$SINGULARITYENV_FOO"}'
        'SINGULARITYENV_X=1'       | 'SINGULARITYENV_X=1'
    }
}
