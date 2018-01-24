/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
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

package nextflow.executor
import java.nio.file.Paths

import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class IgScriptStagingStrategyTest extends Specification {

    def 'should check delegate methods' () {
        given:
        def PATH = Paths.get('/any/path')
        def delegate = new SimpleFileCopyStrategy()
        def strategy = new IgScriptStagingStrategy(delegate: delegate)

        expect:
        strategy.getBeforeStartScript() == null
        strategy.getStageInputFilesScript(['foo.txt': PATH]) == null
        strategy.getUnstageOutputFilesScript(['foo.txt'], PATH) == null
        strategy.resolveForeignFiles(['foo.txt': PATH]) == ['foo.txt': PATH]

        strategy.touchFile(PATH) == delegate.touchFile(PATH)
        strategy.fileStr(PATH) == delegate.fileStr(PATH)
        strategy.copyFile('foo.txt', PATH) == delegate.copyFile('foo.txt', PATH)
        strategy.exitFile(PATH) == delegate.exitFile(PATH)
        strategy.pipeInputFile(PATH) == delegate.pipeInputFile(PATH)
        strategy.getEnvScript([FOO:1],null) == delegate.getEnvScript([FOO:1],null)
    }
}
