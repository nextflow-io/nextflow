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

package test

import java.nio.file.Paths

import com.esotericsoftware.kryo.Kryo
import com.esotericsoftware.kryo.io.Input
import com.esotericsoftware.kryo.io.Output
import nextflow.util.PathSerializer

Kryo kryo = new Kryo();
def clazz = Class.forName("sun.nio.fs.UnixPath")
kryo.register(clazz, new PathSerializer())

def file = File.createTempFile('testkryo',null)
file.deleteOnExit()

def obj = Paths.get('/some/path/hello.txt')
def output = new Output( new FileOutputStream(file) )
kryo.writeClassAndObject(output, obj)
output.close()

def input = new Input( new FileInputStream(file))
def copy = kryo.readClassAndObject(input)
input.close()

println "copy: $copy"

assert copy == Paths.get('/some/path/hello.txt')