/*
 * Copyright 2013-2023, Seqera Labs
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
        new CharliecloudBuilder('busybox')
                .build()
                .runCommand == 'ch-run --unset-env="*" -c "$PWD" --set-env -w -b "$PWD" busybox --'

        new CharliecloudBuilder('busybox')
                .params(runOptions: '-j --no-home')
                .build()
                .runCommand == 'ch-run --unset-env="*" -c "$PWD" --set-env -w -b "$PWD" -j --no-home busybox --'
        
        new CharliecloudBuilder('busybox')
                .params(temp: '/foo')
                .build()
                .runCommand == 'ch-run --unset-env="*" -c "$PWD" --set-env -w -b /foo:/tmp -b "$PWD" busybox --'

        new CharliecloudBuilder('busybox')
                .addEnv('X=1')
                .addEnv(ALPHA:'aaa', BETA: 'bbb')
                .build()
                .runCommand == 'ch-run --unset-env="*" -c "$PWD" --set-env -w --set-env=X=1 --set-env=ALPHA=aaa --set-env=BETA=bbb -b "$PWD" busybox --'

        new CharliecloudBuilder('ubuntu')
                .addMount(path1)
                .build()
                .runCommand == 'ch-run --unset-env="*" -c "$PWD" --set-env -w -b /foo/data/file1 -b "$PWD" ubuntu --'

        new CharliecloudBuilder('ubuntu')
                .addMount(path1)
                .addMount(path2)
                .build()
                .runCommand == 'ch-run --unset-env="*" -c "$PWD" --set-env -w -b /foo/data/file1 -b /bar/data/file2 -b "$PWD" ubuntu --'
    }

    def db_file = Paths.get('/home/db')
    def 'should get run command' () {

        when:
        def cmd = new CharliecloudBuilder('ubuntu').build().getRunCommand()
        then:
        cmd == 'ch-run --unset-env="*" -c "$PWD" --set-env -w -b "$PWD" ubuntu --'

        when:
        cmd = new CharliecloudBuilder('ubuntu').build().getRunCommand('bwa --this --that file.fastq')
        then:
        cmd == 'ch-run --unset-env="*" -c "$PWD" --set-env -w -b "$PWD" ubuntu -- bwa --this --that file.fastq'

        when:
        cmd = new CharliecloudBuilder('ubuntu').params(entry:'/bin/sh').build().getRunCommand('bwa --this --that file.fastq')
        then:
        cmd == 'ch-run --unset-env="*" -c "$PWD" --set-env -w -b "$PWD" ubuntu -- /bin/sh -c "bwa --this --that file.fastq"'

        when:
        cmd = new CharliecloudBuilder('ubuntu').params(entry:'/bin/sh').params(readOnlyInputs: 'true').build().getRunCommand('bwa --this --that file.fastq')
        then:
        cmd == 'ch-run --unset-env="*" -c "$PWD" --set-env -b "$PWD" ubuntu -- /bin/sh -c "bwa --this --that file.fastq"'

        when:
        cmd = new CharliecloudBuilder('ubuntu').params(entry:'/bin/sh').params(readOnlyInputs: 'false').build().getRunCommand('bwa --this --that file.fastq')
        then:
        cmd == 'ch-run --unset-env="*" -c "$PWD" --set-env -w -b "$PWD" ubuntu -- /bin/sh -c "bwa --this --that file.fastq"'

        when:
        cmd = new CharliecloudBuilder('ubuntu')
            .params(entry:'/bin/sh')
            .addMount(db_file)
            .addMount(db_file)
            .params(readOnlyInputs: 'true')
            .build().getRunCommand('bwa --this --that file.fastq')
        then:
        cmd == 'ch-run --unset-env="*" -c "$PWD" --set-env -b /home -b "$PWD" ubuntu -- /bin/sh -c "bwa --this --that file.fastq"'

        when:
        cmd = new CharliecloudBuilder('ubuntu')
            .params(entry:'/bin/sh')
            .addMount(db_file)
            .addMount(db_file)
            .params(readOnlyInputs: 'false')
            .build()
            .getRunCommand('bwa --this --that file.fastq')
        then:
        cmd == 'ch-run --unset-env="*" -c "$PWD" --set-env -w -b /home/db -b "$PWD" ubuntu -- /bin/sh -c "bwa --this --that file.fastq"'
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
