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

import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * @author Patrick Hüther <patrick.huether@gmail.com>
 * @author Laurent Modolo <laurent.modolo@ens-lyon.fr>
 */
class CharliecloudBuilderTest extends Specification {

    def 'should get the ch-run command line' () {

        given:
        def path1 = Paths.get('/foo/data/file1')
        def path2 = Paths.get('/bar/data/file2')

        expect:
        new CharliecloudBuilder('/cacheDir/busybox')
                .build()
                .runCommand == 'ch-run --unset-env="*" -c "$NXF_TASK_WORKDIR" --set-env --write-fake -b "$NXF_TASK_WORKDIR" /cacheDir/busybox --'

        new CharliecloudBuilder('/cacheDir/busybox')
                .params(writeFake: false)
                .build()
                .runCommand == 'ch-run --unset-env="*" -c "$NXF_TASK_WORKDIR" --set-env -w -b "$NXF_TASK_WORKDIR" /cacheDir/busybox --'

        new CharliecloudBuilder('/cacheDir/busybox')
                .params(writeFake: false)
                .params(readOnlyInputs: true)
                .build()
                .runCommand == 'ch-run --unset-env="*" -c "$NXF_TASK_WORKDIR" --set-env -b "$NXF_TASK_WORKDIR" /cacheDir/busybox --'

        new CharliecloudBuilder('/cacheDir/busybox')
                .params(runOptions: '-j --no-home')
                .build()
                .runCommand == 'ch-run --unset-env="*" -c "$NXF_TASK_WORKDIR" --set-env --write-fake -b "$NXF_TASK_WORKDIR" -j --no-home /cacheDir/busybox --'
        
        new CharliecloudBuilder('/cacheDir/busybox')
                .params(temp: '/foo')
                .build()
                .runCommand == 'ch-run --unset-env="*" -c "$NXF_TASK_WORKDIR" --set-env --write-fake -b /foo:/tmp -b "$NXF_TASK_WORKDIR" /cacheDir/busybox --'

        new CharliecloudBuilder('/cacheDir/busybox')
                .addEnv('X=1')
                .addEnv(ALPHA:'aaa', BETA: 'bbb')
                .build()
                .runCommand == 'ch-run --unset-env="*" -c "$NXF_TASK_WORKDIR" --set-env --write-fake --set-env=X=1 --set-env=ALPHA=aaa --set-env=BETA=bbb -b "$NXF_TASK_WORKDIR" /cacheDir/busybox --'

        new CharliecloudBuilder('/cacheDir/ubuntu')
                .addMount(path1)
                .build()
                .runCommand == 'ch-run --unset-env="*" -c "$NXF_TASK_WORKDIR" --set-env --write-fake -b /foo/data/file1 -b "$NXF_TASK_WORKDIR" /cacheDir/ubuntu --'

        new CharliecloudBuilder('/cacheDir/ubuntu')
                .addMount(path1)
                .addMount(path2)
                .build()
                .runCommand == 'ch-run --unset-env="*" -c "$NXF_TASK_WORKDIR" --set-env --write-fake -b /foo/data/file1 -b /bar/data/file2 -b "$NXF_TASK_WORKDIR" /cacheDir/ubuntu --'
    }

    def db_file = Paths.get('/home/db')
    def 'should get run command' () {

        when:
        def cmd = new CharliecloudBuilder('/cacheDir/ubuntu')
            .build()
            .getRunCommand()
        then:
        cmd == 'ch-run --unset-env="*" -c "$NXF_TASK_WORKDIR" --set-env --write-fake -b "$NXF_TASK_WORKDIR" /cacheDir/ubuntu --'

        when:
        cmd = new CharliecloudBuilder('/cacheDir/ubuntu')
            .params(writeFake: 'true')
            .build()
            .getRunCommand()
        then:
        cmd == 'ch-run --unset-env="*" -c "$NXF_TASK_WORKDIR" --set-env --write-fake -b "$NXF_TASK_WORKDIR" /cacheDir/ubuntu --'

        when:
        cmd = new CharliecloudBuilder('/cacheDir/ubuntu')
            .params(writeFake: 'false')
            .build()
            .getRunCommand()
        then:
        cmd == 'ch-run --unset-env="*" -c "$NXF_TASK_WORKDIR" --set-env -w -b "$NXF_TASK_WORKDIR" /cacheDir/ubuntu --'

        when:
        cmd = new CharliecloudBuilder('/cacheDir/ubuntu')
            .params(entry:'/bin/sh')
            .build()
            .getRunCommand('bwa --this --that file.fastq')
        then:
        cmd == 'ch-run --unset-env="*" -c "$NXF_TASK_WORKDIR" --set-env --write-fake -b "$NXF_TASK_WORKDIR" /cacheDir/ubuntu -- /bin/sh -c "bwa --this --that file.fastq"'

        when:
        cmd = new CharliecloudBuilder('/cacheDir/ubuntu')
            .params(entry:'/bin/sh')
            .params(readOnlyInputs: 'true')
            .build()
            .getRunCommand('bwa --this --that file.fastq')
        then:
        cmd == 'ch-run --unset-env="*" -c "$NXF_TASK_WORKDIR" --set-env --write-fake -b "$NXF_TASK_WORKDIR" /cacheDir/ubuntu -- /bin/sh -c "bwa --this --that file.fastq"'

        when:
        cmd = new CharliecloudBuilder('/cacheDir/ubuntu')
            .params(entry:'/bin/sh')
            .params(readOnlyInputs: 'false')
            .build()
            .getRunCommand('bwa --this --that file.fastq')
        then:
        cmd == 'ch-run --unset-env="*" -c "$NXF_TASK_WORKDIR" --set-env --write-fake -b "$NXF_TASK_WORKDIR" /cacheDir/ubuntu -- /bin/sh -c "bwa --this --that file.fastq"'

        when:
        cmd = new CharliecloudBuilder('/cacheDir/ubuntu')
            .params(entry:'/bin/sh')
            .params(readOnlyInputs: 'false')
            .params(writeFake: 'false')
            .build()
            .getRunCommand('bwa --this --that file.fastq')
        then:
        cmd == 'ch-run --unset-env="*" -c "$NXF_TASK_WORKDIR" --set-env -w -b "$NXF_TASK_WORKDIR" /cacheDir/ubuntu -- /bin/sh -c "bwa --this --that file.fastq"'

        when:
        cmd = new CharliecloudBuilder('/cacheDir/ubuntu')
            .params(entry:'/bin/sh')
            .params(readOnlyInputs: 'true')
            .params(writeFake: 'false')
            .addMount(db_file)
            .addMount(db_file)
            .build().getRunCommand('bwa --this --that file.fastq')
        then:
        cmd == 'ch-run --unset-env="*" -c "$NXF_TASK_WORKDIR" --set-env -b /home -b "$NXF_TASK_WORKDIR" /cacheDir/ubuntu -- /bin/sh -c "bwa --this --that file.fastq"'

        when:
        cmd = new CharliecloudBuilder('/cacheDir/ubuntu')
            .params(entry:'/bin/sh')
            .addMount(db_file)
            .addMount(db_file)
            .params(readOnlyInputs: 'false')
            .build()
            .getRunCommand('bwa --this --that file.fastq')
        then:
        cmd == 'ch-run --unset-env="*" -c "$NXF_TASK_WORKDIR" --set-env --write-fake -b /home/db -b "$NXF_TASK_WORKDIR" /cacheDir/ubuntu -- /bin/sh -c "bwa --this --that file.fastq"'
    }

    @Unroll
    def 'test charliecloud env'() {

        given:
        def builder = Spy(CharliecloudBuilder)

        expect:
        builder.makeEnv(ENV).toString() == RESULT

        where:
        ENV                     | RESULT
        'X=1'                   | '--set-env=X=1'
        'BAR'                   | '${BAR:+--set-env=BAR=$BAR}'
        [VAR_X:1, VAR_Y:2]      | '--set-env=VAR_X=1 --set-env=VAR_Y=2'
        [FOO: 'x', BAR: 'y']    | '--set-env=FOO=x --set-env=BAR=y'
    }
}
