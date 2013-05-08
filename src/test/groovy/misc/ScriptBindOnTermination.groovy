package misc
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

import nextflow.Session
/**
 *
 *  @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */


def session = new Session()

def result = session.createProcessor()
        .input(item: [1,2,3])
        .output('result')
        .echo(true)
        .setBindOnTermination(true)
        .setSharedWorkDirectory(true)
        .script {
            """
            echo "item $item" >> result
            """
            }
        .run()


session.createProcessor()
    .echo(true)
    .input(file:result)
    .script { "cat $fInputFile" }
    .run()


session.terminate()