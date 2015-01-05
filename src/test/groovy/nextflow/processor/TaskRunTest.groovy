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

package nextflow.processor

import java.nio.file.Files
import java.nio.file.Paths

import nextflow.file.FileHolder
import nextflow.script.EnvInParam
import nextflow.script.FileInParam
import nextflow.script.FileOutParam
import nextflow.script.FileSharedParam
import nextflow.script.StdInParam
import nextflow.script.StdOutParam
import nextflow.script.TokenVar
import nextflow.script.ValueInParam
import nextflow.script.ValueOutParam
import nextflow.script.ValueSharedParam
import spock.lang.Specification
/**
 *
 *  @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

class TaskRunTest extends Specification {

    def testGetInputsByType() {

        setup:
        def binding = new Binding('x': 1, 'y': 2)
        def task = new TaskRun()
        def list = []

        task.setInput( new StdInParam(binding,list) )
        task.setInput( new FileInParam(binding, list).bind(new TokenVar('x')), 'file1' )
        task.setInput( new FileInParam(binding, list).bind(new TokenVar('y')), 'file2' )
        task.setInput( new EnvInParam(binding, list).bind('z'), 'env' )


        when:
        def files = task.getInputsByType(FileInParam)
        then:
        files.size() == 2

        files.keySet()[0] instanceof FileInParam
        files.keySet()[1] instanceof FileInParam

        files.keySet()[0].name == 'x'
        files.keySet()[1].name == 'y'

        files.values()[0] == 'file1'
        files.values()[1] == 'file2'

    }

    def testGetOutputsByType() {

        setup:
        def binding = new Binding()
        def task = new TaskRun()
        def list = []

        task.setOutput( new FileOutParam(binding, list).bind('x'), 'file1' )
        task.setOutput( new FileOutParam(binding, list).bind('y'), 'file2' )
        task.setOutput( new StdOutParam(binding, list), 'Hello' )


        when:
        def files = task.getOutputsByType(FileOutParam)
        then:
        files.size() == 2

        files.keySet()[0] instanceof FileOutParam
        files.keySet()[1] instanceof FileOutParam

        files.keySet()[0].name == 'x'
        files.keySet()[1].name == 'y'

        files.values()[0] == 'file1'
        files.values()[1] == 'file2'

    }

    def testGetInputFiles() {

        setup:
        def binding = new Binding()
        def task = new TaskRun()
        def list = []

        def x = new ValueInParam(binding, list).bind( new TokenVar('x') )
        def y = new FileInParam(binding, list).bind('y')
        def z = new FileSharedParam(binding, list).bind('z')

        task.setInput(x, 1)
        task.setInput(y, [ new FileHolder(Paths.get('file_y_1')) ])
        task.setInput(z, [ new FileHolder(Paths.get('file_z_1')), new FileHolder(Paths.get('file_z_2')) ])

        expect:
        task.getInputFiles().size() == 2
        task.stagedInputs.size() == 3

    }


    def testGetOutputFilesNames() {

        setup:
        def binding = new Binding()
        def task = new TaskRun()
        def list = []

        when:
        def i1 = new ValueInParam(binding, list).bind( new TokenVar('x') )
        def s1 = new FileSharedParam(binding, list).bind( new TokenVar('y') )
        def o1 = new FileOutParam(binding,list).bind('file_out.alpha')
        def o2 = new ValueOutParam(binding,list).bind( 'x' )
        def o3 = new FileOutParam(binding,list).bind('file_out.beta')

        task.setInput(i1, 'Hello' )
        task.setInput(s1, [ new FileHolder(Paths.get('file_shared.delta')) ])
        task.setOutput(o1)
        task.setOutput(o2)
        task.setOutput(o3)

        then:
        task.getOutputFilesNames() == ['file_out.alpha', 'file_out.beta', 'file_shared.delta']
    }

    def testHasCacheableValues() {

        given:
        def binding = new Binding()
        def list = []

        /*
         * just file as output => no cacheable values
         */
        when:
        def o1 = new FileOutParam(binding,list).bind('file_out.beta')
        def task1 = new TaskRun()
        task1.setOutput(o1)
        then:
        !task1.hasCacheableValues()

        /*
         * when a 'val' is declared as output ==> true
         */
        when:
        def o2 = new ValueOutParam(binding,list).bind( 'x' )
        def task2 = new TaskRun()
        task2.setOutput(o2)
        then:
        task2.hasCacheableValues()

        /*
         * file with parametric name => true
         */
        when:
        def s3 = new FileOutParam(binding, list).bind( new TokenVar('y') )
        def task3 = new TaskRun()
        task3.setOutput(s3)
        then:
        task3.hasCacheableValues()

        /*
         * shared input => true
         */
        when:
        def s4 = new ValueSharedParam(binding, list).bind( 'x' )
        def task4 = new TaskRun()
        task4.setInput(s4)
        then:
        task4.hasCacheableValues()

        when:
        def task5 = new TaskRun( config: new TaskConfig(alpha: 1, beta: 2) )
        then:
        !task5.hasCacheableValues()

        when:
        def task6 = new TaskRun( config: new TaskConfig(alpha: { 'dynamic val' }) )
        then:
        task6.hasCacheableValues()

    }



    def testDumpStdout() {

        setup:
        def file = Files.createTempFile('dumptest',null)
        def task = new TaskRun()

        when:
        task.stdout = '1\n2'
        then:
        task.dumpStdout(5) == ['1','2']

        when:
        task.stdout = """
            1
            2
            3
            4
            5
            6
            7
            8
            9
            """.stripIndent()
        then:
        task.dumpStdout(5) == ['5','6','7','8','9']


        when:
        task = new TaskRun()
        task.stdout = file
        file.text = """
            a
            b
            c
            d
            e
            f
            g
            h
            i
            """.stripIndent()
        then:
        task.dumpStdout(5) == ['e','f','g','h','i']


        when:
        task = new TaskRun()
        task.stdout = Paths.get('/no/file')
        then:
        task.dumpStdout(5) == []

        cleanup:
        file?.delete()

    }

    def testIsSuccess() {

        when:
        def task = new TaskRun(config: [validExitStatus: [0]])
        then:
        task.isSuccess(0) == true
        task.isSuccess('0') == true

        task.isSuccess(1) == false
        task.isSuccess(null) == false
        task.isSuccess('1') == false
        task.isSuccess(Integer.MAX_VALUE) == false


        when:
        task = new TaskRun(config: [validExitStatus: 0])
        then:
        task.isSuccess(0) == true
        task.isSuccess('0') == true

        task.isSuccess(1) == false
        task.isSuccess(null) == false
        task.isSuccess('1') == false
        task.isSuccess(Integer.MAX_VALUE) == false

        when:
        task = new TaskRun(config: [validExitStatus: [0,1]])
        then:
        task.isSuccess(0) == true
        task.isSuccess(1) == true
        task.isSuccess(2) == false

        when:
        task = new TaskRun(config: [validExitStatus: [0]])
        task.exitStatus = 0
        then:
        task.isSuccess() == true

        when:
        task = new TaskRun(config: [validExitStatus: [0]])
        task.exitStatus = 1
        then:
        task.isSuccess() == false

        when:
        task = new TaskRun(config: [validExitStatus: [0]])
        then:
        task.isSuccess() == false
    }


}
