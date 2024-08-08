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
        new CharliecloudBuilder('/cacheDir/img/busybox')
                .params(containerStore: '/containerStorage/')
                .build()
                .runCommand == 'mkdir -p /containerStorage/"${NXF_TASK_WORKDIR: -34}" && ch-convert -i ch-image --storage /cacheDir busybox /containerStorage/"${NXF_TASK_WORKDIR: -34}" && ch-run --unset-env="*" -c "$NXF_TASK_WORKDIR" --set-env -w -b "$NXF_TASK_WORKDIR" /containerStorage/"${NXF_TASK_WORKDIR: -34}" --'

        new CharliecloudBuilder('/cacheDir/img/busybox')
                .params(containerStore: '/containerStorage/')
                .params(runOptions: '-j --no-home')
                .build()
                .runCommand == 'mkdir -p /containerStorage/"${NXF_TASK_WORKDIR: -34}" && ch-convert -i ch-image --storage /cacheDir busybox /containerStorage/"${NXF_TASK_WORKDIR: -34}" && ch-run --unset-env="*" -c "$NXF_TASK_WORKDIR" --set-env -w -b "$NXF_TASK_WORKDIR" -j --no-home /containerStorage/"${NXF_TASK_WORKDIR: -34}" --'
        
        new CharliecloudBuilder('/cacheDir/img/busybox')
                .params(containerStore: '/containerStorage/')
                .params(temp: '/foo')
                .build()
                .runCommand == 'mkdir -p /containerStorage/"${NXF_TASK_WORKDIR: -34}" && ch-convert -i ch-image --storage /cacheDir busybox /containerStorage/"${NXF_TASK_WORKDIR: -34}" && ch-run --unset-env="*" -c "$NXF_TASK_WORKDIR" --set-env -w -b /foo:/tmp -b "$NXF_TASK_WORKDIR" /containerStorage/"${NXF_TASK_WORKDIR: -34}" --'

        new CharliecloudBuilder('/cacheDir/img/busybox')
                .params(containerStore: '/containerStorage/')
                .addEnv('X=1')
                .addEnv(ALPHA:'aaa', BETA: 'bbb')
                .build()
                .runCommand == 'mkdir -p /containerStorage/"${NXF_TASK_WORKDIR: -34}" && ch-convert -i ch-image --storage /cacheDir busybox /containerStorage/"${NXF_TASK_WORKDIR: -34}" && ch-run --unset-env="*" -c "$NXF_TASK_WORKDIR" --set-env -w --set-env=X=1 --set-env=ALPHA=aaa --set-env=BETA=bbb -b "$NXF_TASK_WORKDIR" /containerStorage/"${NXF_TASK_WORKDIR: -34}" --'

        new CharliecloudBuilder('/cacheDir/img/ubuntu')
                .params(containerStore: '/containerStorage/')
                .addMount(path1)
                .build()
                .runCommand == 'mkdir -p /containerStorage/"${NXF_TASK_WORKDIR: -34}" && ch-convert -i ch-image --storage /cacheDir ubuntu /containerStorage/"${NXF_TASK_WORKDIR: -34}" && ch-run --unset-env="*" -c "$NXF_TASK_WORKDIR" --set-env -w -b /foo/data/file1 -b "$NXF_TASK_WORKDIR" /containerStorage/"${NXF_TASK_WORKDIR: -34}" --'

        new CharliecloudBuilder('/cacheDir/img/ubuntu')
                .params(containerStore: '/containerStorage/')
                .addMount(path1)
                .addMount(path2)
                .build()
                .runCommand == 'mkdir -p /containerStorage/"${NXF_TASK_WORKDIR: -34}" && ch-convert -i ch-image --storage /cacheDir ubuntu /containerStorage/"${NXF_TASK_WORKDIR: -34}" && ch-run --unset-env="*" -c "$NXF_TASK_WORKDIR" --set-env -w -b /foo/data/file1 -b /bar/data/file2 -b "$NXF_TASK_WORKDIR" /containerStorage/"${NXF_TASK_WORKDIR: -34}" --'
    }

    def db_file = Paths.get('/home/db')
    def 'should get run command' () {

        when:
        def cmd = new CharliecloudBuilder('/cacheDir/img/ubuntu')
            .params(containerStore: '/containerStorage/')
            .build()
            .getRunCommand()
        then:
        cmd == 'mkdir -p /containerStorage/"${NXF_TASK_WORKDIR: -34}" && ch-convert -i ch-image --storage /cacheDir ubuntu /containerStorage/"${NXF_TASK_WORKDIR: -34}" && ch-run --unset-env="*" -c "$NXF_TASK_WORKDIR" --set-env -w -b "$NXF_TASK_WORKDIR" /containerStorage/"${NXF_TASK_WORKDIR: -34}" --'

        when:
        cmd = new CharliecloudBuilder('/cacheDir/img/ubuntu')
            .params(containerStore: '/containerStorage/')
            .params(useSquash: 'true')
            .build()
            .getRunCommand()
        then:
        cmd == 'mkdir -p /containerStorage/"${NXF_TASK_WORKDIR: -34}" && ch-convert -i ch-image --storage /cacheDir ubuntu /containerStorage/"${NXF_TASK_WORKDIR: -34}"/container.squashfs && ch-run --unset-env="*" -c "$NXF_TASK_WORKDIR" --set-env -w -b "$NXF_TASK_WORKDIR" /containerStorage/"${NXF_TASK_WORKDIR: -34}"/container.squashfs --'
        
        when:
        cmd = new CharliecloudBuilder('/cacheDir/img/ubuntu')
            .params(containerStore: '/containerStorage/')
            .params(writeFake: 'true')
            .build()
            .getRunCommand()
        then:
        cmd == 'ch-run --unset-env="*" -c "$NXF_TASK_WORKDIR" --set-env --write-fake -w -b "$NXF_TASK_WORKDIR" ubuntu --'

        when:
        cmd = new CharliecloudBuilder('/cacheDir/img/ubuntu')
            .params(containerStore: '/containerStorage/')
            .params(entry:'/bin/sh')
            .build()
            .getRunCommand('bwa --this --that file.fastq')
        then:
        cmd == 'mkdir -p /containerStorage/"${NXF_TASK_WORKDIR: -34}" && ch-convert -i ch-image --storage /cacheDir ubuntu /containerStorage/"${NXF_TASK_WORKDIR: -34}" && ch-run --unset-env="*" -c "$NXF_TASK_WORKDIR" --set-env -w -b "$NXF_TASK_WORKDIR" /containerStorage/"${NXF_TASK_WORKDIR: -34}" -- /bin/sh -c "bwa --this --that file.fastq"'

        when:
        cmd = new CharliecloudBuilder('/cacheDir/img/ubuntu')
            .params(containerStore: '/containerStorage/')
            .params(entry:'/bin/sh')
            .params(readOnlyInputs: 'true')
            .build()
            .getRunCommand('bwa --this --that file.fastq')
        then:
        cmd == 'mkdir -p /containerStorage/"${NXF_TASK_WORKDIR: -34}" && ch-convert -i ch-image --storage /cacheDir ubuntu /containerStorage/"${NXF_TASK_WORKDIR: -34}" && ch-run --unset-env="*" -c "$NXF_TASK_WORKDIR" --set-env -b "$NXF_TASK_WORKDIR" /containerStorage/"${NXF_TASK_WORKDIR: -34}" -- /bin/sh -c "bwa --this --that file.fastq"'

        when:
        cmd = new CharliecloudBuilder('/cacheDir/img/ubuntu')
            .params(containerStore: '/containerStorage/')
            .params(entry:'/bin/sh')
            .params(readOnlyInputs: 'false')
            .build()
            .getRunCommand('bwa --this --that file.fastq')
        then:
        cmd == 'mkdir -p /containerStorage/"${NXF_TASK_WORKDIR: -34}" && ch-convert -i ch-image --storage /cacheDir ubuntu /containerStorage/"${NXF_TASK_WORKDIR: -34}" && ch-run --unset-env="*" -c "$NXF_TASK_WORKDIR" --set-env -w -b "$NXF_TASK_WORKDIR" /containerStorage/"${NXF_TASK_WORKDIR: -34}" -- /bin/sh -c "bwa --this --that file.fastq"'

        when:
        cmd = new CharliecloudBuilder('/cacheDir/img/ubuntu')
            .params(containerStore: '/containerStorage/')
            .params(entry:'/bin/sh')
            .addMount(db_file)
            .addMount(db_file)
            .params(readOnlyInputs: 'true')
            .build().getRunCommand('bwa --this --that file.fastq')
        then:
        cmd == 'mkdir -p /containerStorage/"${NXF_TASK_WORKDIR: -34}" && ch-convert -i ch-image --storage /cacheDir ubuntu /containerStorage/"${NXF_TASK_WORKDIR: -34}" && ch-run --unset-env="*" -c "$NXF_TASK_WORKDIR" --set-env -b /home -b "$NXF_TASK_WORKDIR" /containerStorage/"${NXF_TASK_WORKDIR: -34}" -- /bin/sh -c "bwa --this --that file.fastq"'

        when:
        cmd = new CharliecloudBuilder('/cacheDir/img/ubuntu')
            .params(containerStore: '/containerStorage/')
            .params(entry:'/bin/sh')
            .addMount(db_file)
            .addMount(db_file)
            .params(readOnlyInputs: 'false')
            .build()
            .getRunCommand('bwa --this --that file.fastq')
        then:
        cmd == 'mkdir -p /containerStorage/"${NXF_TASK_WORKDIR: -34}" && ch-convert -i ch-image --storage /cacheDir ubuntu /containerStorage/"${NXF_TASK_WORKDIR: -34}" && ch-run --unset-env="*" -c "$NXF_TASK_WORKDIR" --set-env -w -b /home/db -b "$NXF_TASK_WORKDIR" /containerStorage/"${NXF_TASK_WORKDIR: -34}" -- /bin/sh -c "bwa --this --that file.fastq"'
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
