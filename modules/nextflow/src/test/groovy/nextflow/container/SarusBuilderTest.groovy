/*
 * Copyright 2022, Pawsey Supercomputing Research Centre
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
import spock.lang.Specification
import java.nio.file.Paths
/**
 *
 * @author Marco De La Pierre <marco.delapierre@gmail.com>
 */
class SarusrBuilderTest extends Specification {

    def 'test sarus env'() {

        given:
        def builder = new SarusBuilder('x')

        expect:
        builder.makeEnv('X=1').toString() == '-e "X=1"'
        builder.makeEnv([VAR_X:1, VAR_Y: 2]).toString() == '-e "VAR_X=1" -e "VAR_Y=2"'
    }

    def 'should build the sarus run command' () {

        expect:
        new SarusBuilder('busybox')
                .build()
                .@runCommand == 'sarus run --mount=type=bind,source="$PWD",destination="$PWD" -w "$PWD" busybox'

        new SarusBuilder('busybox')
                .params(verbose: true)
                .build()
                .@runCommand == 'sarus --verbose run --mount=type=bind,source="$PWD",destination="$PWD" -w "$PWD" busybox'

        new SarusBuilder('fedora')
                .addEnv([VAR_X:1, VAR_Y:2])
                .addEnv("VAR_Z=3")
                .build()
                .@runCommand == 'sarus run -e "VAR_X=1" -e "VAR_Y=2" -e "VAR_Z=3" --mount=type=bind,source="$PWD",destination="$PWD" -w "$PWD" fedora'

        new SarusBuilder('busybox')
                .params(runOptions: '-x --zeta')
                .build()
                .@runCommand == 'sarus run --mount=type=bind,source="$PWD",destination="$PWD" -w "$PWD" -x --zeta busybox'

        new SarusBuilder('fedora')
                .addEnv([VAR_X:1, VAR_Y:2])
                .addMount(Paths.get('/home/db'))
                .addMount(Paths.get('/home/db'))  // <-- add twice the same to prove that the final string won't contain duplicates
                .build()
                .@runCommand == 'sarus run -e "VAR_X=1" -e "VAR_Y=2" --mount=type=bind,source=/home/db,destination=/home/db --mount=type=bind,source="$PWD",destination="$PWD" -w "$PWD" fedora'

    }

    def 'should get run command line' () {

        when:
        def cli = new SarusBuilder('ubuntu:14').build().getRunCommand()
        then:
        cli ==  '''\
        sarus pull ubuntu:14 1>&2
        sarus run --mount=type=bind,source="$PWD",destination="$PWD" -w "$PWD" ubuntu:14
        '''
        .stripIndent().trim()

        when:
        cli = new SarusBuilder('ubuntu:14').build().getRunCommand('bwa --this --that file.fasta')
        then:
        cli ==  '''\
        sarus pull ubuntu:14 1>&2
        sarus run --mount=type=bind,source="$PWD",destination="$PWD" -w "$PWD" ubuntu:14 bwa --this --that file.fasta
        '''
        .stripIndent().trim()

        when:
        cli = new SarusBuilder('ubuntu:14').params(entry:'/bin/bash').build().getRunCommand('bwa --this --that file.fasta')
        then:
        cli ==  '''\
        sarus pull ubuntu:14 1>&2
        sarus run --mount=type=bind,source="$PWD",destination="$PWD" -w "$PWD" ubuntu:14 /bin/bash -c "bwa --this --that file.fasta"
        '''
        .stripIndent().trim()

    }



}
