/*
 * Copyright 2013-2026, Seqera Labs
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

import nextflow.util.MemoryUnit
import spock.lang.Specification
import spock.lang.Unroll

/**
 * @author Joon-Klaps
 */
class AppleContainerBuilderTest extends Specification {

    def 'test apple mounts'() {

        given:
        def builder = new AppleContainerBuilder('busybox')
        def files  = [Paths.get('/folder/data'), Paths.get('/folder/db'), Paths.get('/folder/db')]
        def real   = [Paths.get('/user/yo/nextflow/bin'), Paths.get('/user/yo/nextflow/work'), Paths.get('/db/pdb/local/data')]
        def quotes = [Paths.get('/folder with blanks/A'), Paths.get('/folder with blanks/B')]

        expect:
        builder.makeVolumes([]).toString()      == '-v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" '
        builder.makeVolumes(files).toString()   == '-v /folder:/folder -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" '
        builder.makeVolumes(real).toString()    == '-v /user/yo/nextflow:/user/yo/nextflow -v /db/pdb/local/data:/db/pdb/local/data -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" '
        builder.makeVolumes(quotes).toString()  == '-v /folder\\ with\\ blanks:/folder\\ with\\ blanks -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" '
        and:
        builder.addMountWorkDir(false).makeVolumes([]).toString()     == ''
        builder.addMountWorkDir(false).makeVolumes(files).toString()  == '-v /folder:/folder '
    }

    @Unroll
    def 'test apple env'() {

        given:
        def builder = new AppleContainerBuilder('busybox')

        expect:
        builder.makeEnv(ENV).toString() == EXPECT

        where:
        ENV                   | EXPECT
        'X=1'                 | '-e "X=1"'
        [VAR_X: 1, VAR_Y: 2] | '-e "VAR_X=1" -e "VAR_Y=2"'
        'BAR'                 | '-e "BAR"'
    }

    def 'test apple create command line'() {

        setup:
        def env     = [FOO: 1, BAR: 'hello world']
        def db_file = Paths.get('/home/db')

        expect:
        new AppleContainerBuilder('fedora')
                .build()
                .runCommand == 'container run -i -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" fedora'

        new AppleContainerBuilder('fedora')
                .addEnv(env)
                .build()
                .runCommand == 'container run -i -e "FOO=1" -e "BAR=hello world" -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" fedora'

        new AppleContainerBuilder('ubuntu', new AppleContainerConfig(temp: '/hola'))
                .build()
                .runCommand == 'container run -i -v /hola:/tmp -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" ubuntu'

        new AppleContainerBuilder('busybox')
                .params(entry: '/bin/bash')
                .build()
                .runCommand == 'container run -i -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" --entrypoint /bin/bash busybox'

        new AppleContainerBuilder('busybox', new AppleContainerConfig(runOptions: '-x --zeta'))
                .build()
                .runCommand == 'container run -i -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" -x --zeta busybox'

        new AppleContainerBuilder('busybox')
                .setName('hola')
                .build()
                .runCommand == 'container run -i -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" --name hola busybox'

        new AppleContainerBuilder('busybox', new AppleContainerConfig(engineOptions: '--some-engine-opt'))
                .build()
                .runCommand == 'container --some-engine-opt run -i -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" busybox'

        new AppleContainerBuilder('fedora')
                .addEnv(env)
                .addMount(db_file)
                .addMount(db_file)  // add twice to prove deduplication
                .build()
                .runCommand == 'container run -i -e "FOO=1" -e "BAR=hello world" -v /home/db:/home/db -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" fedora'

        new AppleContainerBuilder('fedora')
                .params(readOnlyInputs: true)
                .addMount(db_file)
                .build()
                .runCommand == 'container run -i -v /home/db:/home/db:ro -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" fedora'
    }

    def 'test add mount'() {

        when:
        def apple = new AppleContainerBuilder('fedora')
        apple.addMount(Paths.get('hello'))
        then:
        apple.mounts.size() == 1
        apple.mounts.contains(Paths.get('hello'))

        when:
        apple.addMount(null)
        then:
        apple.mounts.size() == 1
    }

    def 'test get commands'() {

        when:
        def apple = new AppleContainerBuilder('busybox').setName('c1').build()
        then:
        apple.runCommand    == 'container run -i -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" --name c1 busybox'
        apple.removeCommand == 'container rm c1'
        apple.killCommand   == 'container stop c1'

        when:
        def config = new AppleContainerConfig(remove: true)
        apple = new AppleContainerBuilder('busybox', config).setName('c3').build()
        then:
        apple.runCommand    == 'container run -i -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" --name c3 busybox'
        apple.removeCommand == 'container rm c3'
        apple.killCommand   == 'container stop c3'

        when:
        apple = new AppleContainerBuilder('busybox').setName('c4').params(kill: 'SIGKILL').build()
        then:
        apple.killCommand == 'container kill --signal SIGKILL c4'

        when:
        config = new AppleContainerConfig(remove: false)
        apple = new AppleContainerBuilder('busybox', config).setName('c5').params(kill: false).build()
        then:
        apple.killCommand   == null
        apple.removeCommand == null

        when:
        // no name set — kill and remove commands must be null regardless
        apple = new AppleContainerBuilder('busybox').build()
        then:
        apple.killCommand   == null
        apple.removeCommand == null
    }

    def 'should get run command line'() {

        when:
        def cli = new AppleContainerBuilder('ubuntu:14').build().getRunCommand()
        then:
        cli == 'container run -i -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" ubuntu:14'

        when:
        cli = new AppleContainerBuilder('ubuntu:14').build().getRunCommand(null)
        then:
        cli == 'container run -i -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" ubuntu:14'

        when:
        cli = new AppleContainerBuilder('ubuntu:14').build().getRunCommand('bwa --this --that file.fasta')
        then:
        cli == 'container run -i -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" ubuntu:14 bwa --this --that file.fasta'

        when:
        cli = new AppleContainerBuilder('ubuntu:14').params(entry: '/bin/bash').build().getRunCommand('bwa --this --that file.fasta')
        then:
        cli == 'container run -i -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" --entrypoint /bin/bash ubuntu:14 -c "bwa --this --that file.fasta"'
    }

    def 'test cpus and memory'() {

        expect:
        // --cpus N (not --cpu-shares as in Docker/Podman)
        new AppleContainerBuilder('fedora')
                .setCpus(3)
                .build()
                .runCommand == 'container run -i -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" --cpus 3 fedora'

        // memory uses uppercase suffix as required by the Apple container CLI
        new AppleContainerBuilder('fedora')
                .setMemory(new MemoryUnit('100m'))
                .build()
                .runCommand == 'container run -i -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" --memory 100M fedora'

        new AppleContainerBuilder('fedora')
                .setCpus(1)
                .setMemory(new MemoryUnit('400m'))
                .build()
                .runCommand == 'container run -i -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" --cpus 1 --memory 400M fedora'
    }

    def 'test platform and rosetta'() {

        expect:
        // no platform flag by default
        new AppleContainerBuilder('fedora')
                .build()
                .runCommand == 'container run -i -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" fedora'

        // native arm64 — pass --platform explicitly
        new AppleContainerBuilder('fedora')
                .setPlatform('linux/arm64')
                .build()
                .runCommand == 'container run -i -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" --platform linux/arm64 fedora'

        // amd64 image: --platform selects the image variant, --rosetta enables x86_64 execution
        new AppleContainerBuilder('fedora')
                .setPlatform('linux/amd64')
                .build()
                .runCommand == 'container run -i -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" --platform linux/amd64 --rosetta fedora'

        new AppleContainerBuilder('fedora')
                .setPlatform('amd64')
                .build()
                .runCommand == 'container run -i -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" --platform amd64 --rosetta fedora'
    }

    def 'test registry prefix'() {

        expect:
        new AppleContainerBuilder('ubuntu', new AppleContainerConfig(registry: 'my.registry.io/'))
                .build()
                .runCommand == 'container run -i -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" my.registry.io/ubuntu'
    }

    def 'test factory creates apple builder'() {

        when:
        def config  = new AppleContainerConfig(enabled: true)
        def builder = ContainerBuilder.create(config, 'ubuntu:22.04')

        then:
        builder instanceof AppleContainerBuilder
        builder.image == 'ubuntu:22.04'
    }
}
