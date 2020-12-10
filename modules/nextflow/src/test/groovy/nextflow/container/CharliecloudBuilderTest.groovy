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
 * @author Patrick Hüther <patrick.huether@gmail.com>
 */
class CharliecloudBuilderTest extends Specification {

    def 'should get the ch-run command line' () {

        given:
        def path1 = Paths.get('/foo/data/file1')
        def path2 = Paths.get('/bar/data/file2')

        expect:
        new CharliecloudBuilder('busybox')
                .build()
                .runCommand == 'ch-run --no-home -w busybox -- bash -c "mkdir -p "$PWD"";ch-run --no-home --unset-env="*" -b "$PWD":"$PWD" -c "$PWD" busybox --'

        new CharliecloudBuilder('busybox')
                .params(runOptions: '-j -w')
                .build()
                .runCommand == 'ch-run --no-home -w busybox -- bash -c "mkdir -p "$PWD"";ch-run --no-home --unset-env="*" -j -w -b "$PWD":"$PWD" -c "$PWD" busybox --'

        new CharliecloudBuilder('busybox')
                .addEnv('X=1')
                .addEnv(ALPHA:'aaa', BETA: 'bbb')
                .build()
                .runCommand == 'ch-run --no-home -w busybox -- bash -c "mkdir -p "$PWD"";ch-run --no-home --unset-env="*" --set-env=<( echo "X=1" ) --set-env=<( echo "ALPHA="aaa"" ) --set-env=<( echo "BETA="bbb"" ) -b "$PWD":"$PWD" -c "$PWD" busybox --'

        new CharliecloudBuilder('ubuntu')
                .addMount(path1)
                .build()
                .runCommand == 'ch-run --no-home -w ubuntu -- bash -c "mkdir -p /foo/data/file1 "$PWD"";ch-run --no-home --unset-env="*" -b /foo/data/file1:/foo/data/file1 -b "$PWD":"$PWD" -c "$PWD" ubuntu --'

        new CharliecloudBuilder('ubuntu')
                .addMount(path1)
                .addMount(path2)
                .build()
                .runCommand == 'ch-run --no-home -w ubuntu -- bash -c "mkdir -p /foo/data/file1 /bar/data/file2 "$PWD"";ch-run --no-home --unset-env="*" -b /foo/data/file1:/foo/data/file1 -b /bar/data/file2:/bar/data/file2 -b "$PWD":"$PWD" -c "$PWD" ubuntu --'

        new CharliecloudBuilder('ubuntu')
                .addMount(path1)
                .addMount(path2)
                .params(readOnlyInputs: true)
                .build()
                .runCommand == 'ch-run --no-home -w ubuntu -- bash -c "mkdir -p /foo/data/file1 /bar/data/file2 "$PWD"";ch-run --no-home --unset-env="*" -b /foo/data/file1:/foo/data/file1 -b /bar/data/file2:/bar/data/file2 -b "$PWD":"$PWD" -c "$PWD" ubuntu --'
    }

    def 'should get run command' () {

        when:
        def cmd = new CharliecloudBuilder('ubuntu').build().getRunCommand()
        then:
        cmd == 'ch-run --no-home -w ubuntu -- bash -c "mkdir -p "$PWD"";ch-run --no-home --unset-env="*" -b "$PWD":"$PWD" -c "$PWD" ubuntu --'

        when:
        cmd = new CharliecloudBuilder('ubuntu').build().getRunCommand('bwa --this --that file.fastq')
        then:
        cmd == 'ch-run --no-home -w ubuntu -- bash -c "mkdir -p "$PWD"";ch-run --no-home --unset-env="*" -b "$PWD":"$PWD" -c "$PWD" ubuntu -- bwa --this --that file.fastq'

        when:
        cmd = new CharliecloudBuilder('ubuntu').params(entry:'/bin/sh').build().getRunCommand('bwa --this --that file.fastq')
        then:
        cmd == 'ch-run --no-home -w ubuntu -- bash -c "mkdir -p "$PWD"";ch-run --no-home --unset-env="*" -b "$PWD":"$PWD" -c "$PWD" ubuntu -- /bin/sh -c "bwa --this --that file.fastq"'
    }

    @Unroll
    def 'test charliecloud env'() {

        given:
        def builder = Spy(CharliecloudBuilder)

        expect:
        builder.makeEnv(ENV).toString() == RESULT

        where:
        ENV                     | RESULT
        'X=1'                   | '--set-env=<( echo "X=1" )'
        'BAR'                   | '${BAR:+--set-env=<( echo "BAR="$BAR"" )}'
        [VAR_X:1, VAR_Y:2]      | '--set-env=<( echo "VAR_X="1"" ) --set-env=<( echo "VAR_Y="2"" )'
        [FOO: 'x', BAR: 'y']    | '--set-env=<( echo "FOO="x"" ) --set-env=<( echo "BAR="y"" )'
    }
}
