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

import groovyx.gpars.dataflow.Dataflow
import groovyx.gpars.dataflow.DataflowVariable

/**
 *
 *  @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */


def a = new DataflowVariable()
def b = new DataflowVariable()
def c = new DataflowVariable()
def d = new DataflowVariable()


Dataflow.operator(inputs: [a, b, c], outputs: [d], maxForks: 3) {x, y, z ->
    println "before bind"
    bindOutput 0, x + y + z
    println "after bind"
}


a << 1
b << 2
c << 3

println "before print"
println "=> ${d.val}"

println("Done")