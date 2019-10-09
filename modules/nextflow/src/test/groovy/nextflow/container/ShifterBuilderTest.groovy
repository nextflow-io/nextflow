/*
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
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ShifterBuilderTest extends Specification {

    def 'test shifter env'() {

        given:
        def builder = new ShifterBuilder('x')

        expect:
        builder.makeEnv('X=1').toString() == 'X=1'
        builder.makeEnv([VAR_X:1, VAR_Y: 2]).toString() == 'VAR_X="1" VAR_Y="2"'
    }

    def 'should build the shifter run command' () {

        expect:
        new ShifterBuilder('busybox')
                .build()
                .@runCommand == 'shifter --image busybox'

        new ShifterBuilder('busybox')
                .params(verbose: true)
                .build()
                .@runCommand == 'shifter --verbose --image busybox'

        new ShifterBuilder('fedora')
                .addEnv([VAR_X:1, VAR_Y:2])
                .addEnv("VAR_Z=3")
                .build()
                .@runCommand == 'VAR_X="1" VAR_Y="2" VAR_Z=3 shifter --image fedora'

    }

    def 'should get run command line' () {

        when:
        def cli = new ShifterBuilder('ubuntu:14').build().getRunCommand()
        then:
        cli ==  '''\
        shifterimg pull ubuntu:14
        shifterimg lookup ubuntu:14
        while ! shifterimg lookup ubuntu:14; do
            sleep 5
            STATUS=$(shifterimg -v pull ubuntu:14 | tail -n2 | head -n1 | awk '{print $6}')
            [[ $STATUS == "FAILURE" || -z $STATUS ]] && echo "Shifter failed to pull image 'ubuntu:14'" >&2  && exit 1
        done
        shifter --image ubuntu:14
        '''
        .stripIndent().trim()

        when:
        cli = new ShifterBuilder('ubuntu:14').build().getRunCommand('bwa --this --that file.fasta')
        then:
        cli ==  '''\
        shifterimg pull ubuntu:14
        shifterimg lookup ubuntu:14
        while ! shifterimg lookup ubuntu:14; do
            sleep 5
            STATUS=$(shifterimg -v pull ubuntu:14 | tail -n2 | head -n1 | awk '{print $6}')
            [[ $STATUS == "FAILURE" || -z $STATUS ]] && echo "Shifter failed to pull image 'ubuntu:14'" >&2  && exit 1
        done
        shifter --image ubuntu:14 bwa --this --that file.fasta
        '''
        .stripIndent().trim()

        when:
        cli = new ShifterBuilder('ubuntu:14').params(entry:'/bin/bash').build().getRunCommand('bwa --this --that file.fasta')
        then:
        cli ==  '''\
        shifterimg pull ubuntu:14
        shifterimg lookup ubuntu:14
        while ! shifterimg lookup ubuntu:14; do
            sleep 5
            STATUS=$(shifterimg -v pull ubuntu:14 | tail -n2 | head -n1 | awk '{print $6}')
            [[ $STATUS == "FAILURE" || -z $STATUS ]] && echo "Shifter failed to pull image 'ubuntu:14'" >&2  && exit 1
        done
        shifter --image ubuntu:14 /bin/bash -c "bwa --this --that file.fasta"
        '''
        .stripIndent().trim()

    }



}
