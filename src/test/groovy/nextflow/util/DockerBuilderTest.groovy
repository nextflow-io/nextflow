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



    def testDockerMounts() {

        setup:
        def files =  [Paths.get('/folder/data'),  Paths.get('/folder/db'), Paths.get('/folder/db') ]
        def real = [ Paths.get('/user/yo/nextflow/bin'), Paths.get('/user/yo/nextflow/work'), Paths.get('/db/pdb/local/data') ]

        expect:
        DockerBuilder.makeVolumes([]).toString() == '-v $PWD:$PWD'
        DockerBuilder.makeVolumes(files).toString() == '-v /folder:/folder -v $PWD:$PWD'
        DockerBuilder.makeVolumes(real).toString()  == '-v /user/yo/nextflow:/user/yo/nextflow -v /db/pdb/local/data:/db/pdb/local/data -v $PWD:$PWD'
    }


    def testDockerEnv() {

        expect:
        DockerBuilder.makeEnv('X=1').toString() == '-e "X=1"'
        DockerBuilder.makeEnv([VAR_X:1, VAR_Y: 2]).toString() == '-e "VAR_X=1" -e "VAR_Y=2"'
        DockerBuilder.makeEnv( Paths.get('/some/file.env') ).toString() == '-e "BASH_ENV=file.env"'
        DockerBuilder.makeEnv( new File('/some/file.env') ).toString() == '-e "BASH_ENV=file.env"'
    }

    def testDockerRunCommandLine() {

        setup:
        def envFile = Paths.get('env-file')
        def db_file = Paths.get('/home/db')

        expect:
        new DockerBuilder('fedora').build().toString() == 'docker run --rm -v $PWD:$PWD -w $PWD fedora'
        new DockerBuilder('fedora').addEnv(envFile).build().toString() == 'docker run --rm -e "BASH_ENV=env-file" -v $PWD:$PWD -w $PWD fedora'
        new DockerBuilder('ubuntu').params(temp:'/hola').build().toString() == 'docker run --rm -v /hola:/tmp -v $PWD:$PWD -w $PWD ubuntu'
        new DockerBuilder('ubuntu').params(rm:false).build().toString() == 'docker run -v $PWD:$PWD -w $PWD ubuntu'
        new DockerBuilder('busybox').params(sudo: true).build().toString() == 'sudo docker run --rm -v $PWD:$PWD -w $PWD busybox'
        new DockerBuilder('busybox').params(options: '-x --zeta').build().toString() == 'docker run --rm -v $PWD:$PWD -w $PWD -x --zeta busybox'
        new DockerBuilder('busybox').params(userEmulation:true).build().toString() == 'docker run --rm -u $(id -u) -e "HOME=${HOME}" -v /etc/passwd:/etc/passwd:ro -v /etc/shadow:/etc/shadow:ro -v /etc/group:/etc/group:ro -v $HOME:$HOME -v $PWD:$PWD -w $PWD busybox'
        new DockerBuilder('fedora')
                .addEnv(envFile)
                .addMount(db_file)
                .addMount(db_file) // <-- add twice the same to prove that the final string won't contain duplicates
                .build()
                .toString() == 'docker run --rm -e "BASH_ENV=env-file" -v /home/db:/home/db -v $PWD:$PWD -w $PWD fedora'

    }

    def testAddMount() {

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


}
