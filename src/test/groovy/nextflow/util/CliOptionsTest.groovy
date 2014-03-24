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

import nextflow.script.CliOptions
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CliOptionsTest extends Specification {


    def testNormalizeCmdline () {

        expect:
        CliOptions.normalizeArgs('a','-bb','-ccc','dddd') == ['a','-bb','-ccc','dddd']
        CliOptions.normalizeArgs('a','-bb','-ccc','-resume', 'last') == ['a','-bb','-ccc','-resume','last']
        CliOptions.normalizeArgs('a','-bb','-ccc','-resume') == ['a','-bb','-ccc','-resume','last']
        CliOptions.normalizeArgs('a','-bb','-ccc','-resume','1d2c942a-345d-420b-b7c7-18d90afc6c33', 'zzz') == ['a','-bb','-ccc','-resume','1d2c942a-345d-420b-b7c7-18d90afc6c33', 'zzz']

        CliOptions.normalizeArgs('x','-test') == ['x','-test','%all']
        CliOptions.normalizeArgs('x','-test','alpha') == ['x','-test','alpha']
        CliOptions.normalizeArgs('x','-test','-other') == ['x','-test','%all','-other']

        CliOptions.normalizeArgs('--alpha=1') == ['--alpha=1']
        CliOptions.normalizeArgs('--alpha','1') == ['--alpha=1']
        CliOptions.normalizeArgs('-x', '1', 'script.nf', '--long', 'v1', '--more', 'v2', '--flag') == ['-x','1','script.nf','--long=v1','--more=v2','--flag=true']

        CliOptions.normalizeArgs('-x', '1', '-process.alpha','2', '3') == ['-x', '1', '-process.alpha=2', '3']
        CliOptions.normalizeArgs('-x', '1', '-process.echo') == ['-x', '1', '-process.echo=true']

        CliOptions.normalizeArgs('-x', '1', '-daemon.alpha','2', '3') == ['-x', '1', '-daemon.alpha=2', '3']
        CliOptions.normalizeArgs('-x', '1', '-daemon.echo') == ['-x', '1', '-daemon.echo=true']

        CliOptions.normalizeArgs('-x', '1', '-executor.alpha','2', '3') == ['-x', '1', '-executor.alpha=2', '3']
        CliOptions.normalizeArgs('-x', '1', '-executor.echo') == ['-x', '1', '-executor.echo=true']

        CliOptions.normalizeArgs('-x', '1', '-that.alpha','2', '3') == ['-x', '1', '-that.alpha','2', '3']
    }


}
