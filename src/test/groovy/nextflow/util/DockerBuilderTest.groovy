/*
 * Copyright (c) 2013-2014, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2014, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
 */

package nextflow.util
import java.nio.file.Paths

import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DockerBuilderTest extends Specification {


    def 'test docker mounts'() {

        setup:
        def files =  [Paths.get('/folder/data'),  Paths.get('/folder/db'), Paths.get('/folder/db') ]
        def real = [ Paths.get('/user/yo/nextflow/bin'), Paths.get('/user/yo/nextflow/work'), Paths.get('/db/pdb/local/data') ]

        expect:
        DockerBuilder.makeVolumes([]).toString() == '-v $PWD:$PWD'
        DockerBuilder.makeVolumes(files).toString() == '-v /folder:/folder -v $PWD:$PWD'
        DockerBuilder.makeVolumes(real).toString()  == '-v /user/yo/nextflow:/user/yo/nextflow -v /db/pdb/local/data:/db/pdb/local/data -v $PWD:$PWD'
    }


    def 'test docker env'() {

        expect:
        DockerBuilder.makeEnv('X=1').toString() == '-e "X=1"'
        DockerBuilder.makeEnv([VAR_X:1, VAR_Y: 2]).toString() == '-e "VAR_X=1" -e "VAR_Y=2"'
        DockerBuilder.makeEnv( Paths.get('/some/file.env') ).toString() == '-e "BASH_ENV=file.env"'
        DockerBuilder.makeEnv( new File('/some/file.env') ).toString() == '-e "BASH_ENV=file.env"'
    }

    def 'test docker run command line'() {

        setup:
        def envFile = Paths.get('env-file')
        def db_file = Paths.get('/home/db')

        expect:
        new DockerBuilder('fedora').build() == 'docker run -i -v $PWD:$PWD -w $PWD fedora'
        new DockerBuilder('fedora').addEnv(envFile).build() == 'docker run -i -e "BASH_ENV=env-file" -v $PWD:$PWD -w $PWD fedora'
        new DockerBuilder('ubuntu').params(temp:'/hola').build() == 'docker run -i -v /hola:/tmp -v $PWD:$PWD -w $PWD ubuntu'
        new DockerBuilder('busybox').params(sudo: true).build() == 'sudo docker run -i -v $PWD:$PWD -w $PWD busybox'
        new DockerBuilder('busybox').params(entry: '/bin/bash').build() == 'docker run -i -v $PWD:$PWD -w $PWD --entrypoint /bin/bash busybox'
        new DockerBuilder('busybox').params(runOptions: '-x --zeta').build() == 'docker run -i -v $PWD:$PWD -w $PWD -x --zeta busybox'
        new DockerBuilder('busybox').params(userEmulation:true).build() == 'docker run -i -u $(id -u) -e "HOME=${HOME}" -v /etc/passwd:/etc/passwd:ro -v /etc/shadow:/etc/shadow:ro -v /etc/group:/etc/group:ro -v $HOME:$HOME -v $PWD:$PWD -w $PWD busybox'
        new DockerBuilder('busybox').setName('hola').build() == 'docker run -i -v $PWD:$PWD -w $PWD --name hola busybox'

        new DockerBuilder('fedora')
                .addEnv(envFile)
                .addMount(db_file)
                .addMount(db_file) // <-- add twice the same to prove that the final string won't contain duplicates
                .build() == 'docker run -i -e "BASH_ENV=env-file" -v /home/db:/home/db -v $PWD:$PWD -w $PWD fedora'

    }

    def 'test memory and cpuset' () {

        expect:
        new DockerBuilder('fedora').setCpus('1,2').build() == 'docker run -i --cpuset 1,2 -v $PWD:$PWD -w $PWD fedora'
        new DockerBuilder('fedora').setMemory('10g').build() == 'docker run -i --memory 10g -v $PWD:$PWD -w $PWD fedora'
        new DockerBuilder('fedora').setMemory(new MemoryUnit('100M')).build() == 'docker run -i --memory 100m -v $PWD:$PWD -w $PWD fedora'
        new DockerBuilder('fedora').setCpus('1-3').setMemory(new MemoryUnit('100M')).build() == 'docker run -i --cpuset 1-3 --memory 100m -v $PWD:$PWD -w $PWD fedora'

    }

    def 'test add mount'() {

        when:
        def docker = new DockerBuilder('fedora')
        docker.addMount(Paths.get('hello'))
        then:
        docker.mounts.size() == 1
        docker.mounts.contains(Paths.get('hello'))

        when:
        docker.addMount(null)
        then:
        docker.mounts.size() == 1

    }

    def 'test get commands'() {

        when:
        def docker =  new DockerBuilder('busybox').setName('c1')
        then:
        docker.build() == 'docker run -i -v $PWD:$PWD -w $PWD --name c1 busybox'
        docker.runCommand == 'docker run -i -v $PWD:$PWD -w $PWD --name c1 busybox'
        docker.removeCommand == null
        docker.killCommand == 'docker kill c1'


        when:
        docker =  new DockerBuilder('busybox').setName('c2').params(sudo: true, remove: true)
        then:
        docker.build() == 'sudo docker run -i -v $PWD:$PWD -w $PWD --name c2 busybox'
        docker.runCommand == 'sudo docker run -i -v $PWD:$PWD -w $PWD --name c2 busybox'
        docker.removeCommand == 'sudo docker rm c2'
        docker.killCommand == 'sudo docker kill c2'


        when:
        docker =  new DockerBuilder('busybox').setName('c3').params(remove: true)
        then:
        docker.build() == 'docker run -i -v $PWD:$PWD -w $PWD --name c3 busybox'
        docker.runCommand == 'docker run -i -v $PWD:$PWD -w $PWD --name c3 busybox'
        docker.removeCommand == 'docker rm c3'
        docker.killCommand == 'docker kill c3'

    }

    def 'test is absolute image name' () {

        expect:
        !DockerBuilder.isAbsoluteDockerName('hello')
        !DockerBuilder.isAbsoluteDockerName('image/name')
        DockerBuilder.isAbsoluteDockerName('registry:5000/image/name')
        DockerBuilder.isAbsoluteDockerName('d.reg/image/name')
        DockerBuilder.isAbsoluteDockerName('d.reg/image')

    }

    def 'test normalize docker image name' () {

        expect:
        DockerBuilder.normalizeDockerImageName(image, [registry: registry]) == expected

        where:
        image                       | registry  | expected
        null                        | null      | null
        null                        | 'd.reg'   | null
        'hello'                     | null      | 'hello'
        'cbcrg/hello'               | null      | 'cbcrg/hello'
        'cbcrg/hello'               | 'd.reg'   | 'd.reg/cbcrg/hello'
        'cbcrg/hello'               | 'd.reg/'  | 'd.reg/cbcrg/hello'
        'registry:5000/cbcrg/hello' | 'd.reg'   | 'registry:5000/cbcrg/hello'

    }


}
