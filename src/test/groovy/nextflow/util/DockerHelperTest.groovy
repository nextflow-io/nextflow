/*
 * Copyright (c) 2012, the authors.
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
class DockerHelperTest extends Specification {



    def testDockerMounts() {

        setup:
        def files =  [Paths.get('/folder/data'),  Paths.get('/folder/db'), Paths.get('/folder/db') ]
        def real = [ Paths.get('/user/yo/nextflow/bin'), Paths.get('/user/yo/nextflow/work'), Paths.get('/db/pdb/local/data') ]

        expect:
        DockerHelper.getVolumes([]).toString() == '-v $PWD:$PWD'
        DockerHelper.getVolumes(files).toString() == '-v /folder:/folder -v $PWD:$PWD'
        DockerHelper.getVolumes(real).toString()  == '-v /user/yo/nextflow:/user/yo/nextflow -v /db/pdb/local/data:/db/pdb/local/data -v $PWD:$PWD'
    }


    def testDockerEnv() {

        expect:
        DockerHelper.getEnv('X=1').toString() == '-e "X=1"'
        DockerHelper.getEnv([VAR_X:1, VAR_Y: 2]).toString() == '-e "VAR_X=1" -e "VAR_Y=2"'
        DockerHelper.getEnv( Paths.get('/some/file.env') ).toString() == '-e "BASH_ENV=file.env"'
        DockerHelper.getEnv( new File('/some/file.env') ).toString() == '-e "BASH_ENV=file.env"'
    }

    def testDockerRunCommandLine() {

        setup:
        def envFile = Paths.get('env-file')
        def files =  [Paths.get('/home/db'), Paths.get('/home/db') ]

        expect:
        DockerHelper.getRun('fedora', [], null ).toString() == 'docker run --rm -v $PWD:$PWD -w $PWD fedora'
        DockerHelper.getRun('fedora', [], envFile ).toString() == 'docker run --rm -e "BASH_ENV=env-file" -v $PWD:$PWD -w $PWD fedora'
        DockerHelper.getRun('fedora', files, envFile).toString() == 'docker run --rm -e "BASH_ENV=env-file" -v /home/db:/home/db -v $PWD:$PWD -w $PWD fedora'

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
