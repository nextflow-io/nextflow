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

package nextflow.processor
import nextflow.script.FileInParam
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class MergeTaskProcessorTest extends Specification {

    def testProvider () {

        setup:
        def binding = new Binding()
        def files = [:]
        def key1 = new FileInParam(binding, []).bind('file1')
        def key2 = new FileInParam(binding, []).bind('file_')
        def val1 = [ FileHolder.get('xxx', 'file1') ]
        def val2 =  [ FileHolder.get('yyy', 'file2'), FileHolder.get('zzz', 'file3') ]
        def val3 =  'just a value'
        files[key1] = val1
        files[key2] = val2


        def proc = [:] as MergeTaskProcessor
        proc.filesCollector = files

        expect:
        proc.stagedProvider().size() == 2
        proc.stagedProvider().get(key1) == val1
        proc.stagedProvider().get(key2) == val2


    }

}
