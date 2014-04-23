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
        DockerBuilder.makeVolumes([]).toString() == '-v ${NXF_SCRATCH:-$(mktemp -d)}:/tmp -v $PWD:$PWD'
        DockerBuilder.makeVolumes(files).toString() == '-v ${NXF_SCRATCH:-$(mktemp -d)}:/tmp -v /folder:/folder -v $PWD:$PWD'
        DockerBuilder.makeVolumes(real).toString()  == '-v ${NXF_SCRATCH:-$(mktemp -d)}:/tmp -v /user/yo/nextflow:/user/yo/nextflow -v /db/pdb/local/data:/db/pdb/local/data -v $PWD:$PWD'
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
        def files =  [Paths.get('/home/db'), Paths.get('/home/db') ]

        expect:
        new DockerBuilder('fedora').build().toString() == 'docker run --rm -u $(id -u) -v ${NXF_SCRATCH:-$(mktemp -d)}:/tmp -v $PWD:$PWD -w $PWD fedora'
        new DockerBuilder('fedora').addEnv(envFile).build().toString() == 'docker run --rm -u $(id -u) -e "BASH_ENV=env-file" -v ${NXF_SCRATCH:-$(mktemp -d)}:/tmp -v $PWD:$PWD -w $PWD fedora'
        new DockerBuilder('fedora').addEnv(envFile).addMount(files).build().toString() == 'docker run --rm -u $(id -u) -e "BASH_ENV=env-file" -v ${NXF_SCRATCH:-$(mktemp -d)}:/tmp -v /home/db:/home/db -v $PWD:$PWD -w $PWD fedora'

    }


//
//    def testDockerRunCommandLineWithFiles() {
//
//        when:
//        def fileHolders = new HashMap<?,List<FileHolder>>()
//        fileHolders[ new FileInParam(null,[]) ] = [FileHolder.get('/home/data/sequences', 'file.txt')]
//        fileHolders[ new FileInParam(null,[]) ] = [FileHolder.get('/home/data/file1','seq_1.fa'), FileHolder.get('/home/data/file2','seq_2.fa'), FileHolder.get('/home/data/file3','seq_3.fa') ]
//
//        then:
//        DockerHelper.getRun('ubuntu', fileHolders, null ).toString() == 'docker run --rm -v /home/data:/home/data -v $PWD:$PWD -w $PWD ubuntu'
//
//    }


}
