
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

/**
 *
 *  @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */


def x=1;
def y=2
def z=3

def out1 = 1

def task( Closure<String> closure ) {


    closure.delegate = new MyDelegate()
    def result = closure.call()

    println result
}

class MyDelegate implements GroovyObject {

    Object invokeMethod(String name, Object args) {
        println "method: $name ($args)"
        return 'x'
    }

    Object getProperty(String name) {
        return 'y'
    }

}


task {

    stdin(channelIn)
    stdout(channelOut)

    """
    do this ${input(channelIn)}
    do that >> ${outfile(file:'dssd' )}
               ${out(channel3:'resutl_*', store:file)}
    do that ${shell(variable)}
    """

}


