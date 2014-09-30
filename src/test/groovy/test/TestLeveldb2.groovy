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

package test
import static org.iq80.leveldb.impl.Iq80DBFactory.factory

import java.nio.ByteBuffer

import groovy.transform.CompileStatic
import nextflow.util.KryoHelper
import org.iq80.leveldb.Options

@CompileStatic
class Entry implements Serializable {
    int pos
    String value

    Entry() {}
}

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

buffer = ByteBuffer.allocate(4)

def byte[] bytes(int value) {
    buffer.clear()
    buffer.putInt(value).array();
}

def int asInt(byte[] bytes) {
    ByteBuffer.wrap(bytes).getInt()
}

int count = 0
final EMPTY = [] as byte[]

Options options = new Options();
options.createIfMissing(true);
data = factory.open(new File("/Users/pditommaso/Downloads/db/data"), options);
index = factory.open(new File("/Users/pditommaso/Downloads/db/index"), options);

println "Begin"
def now = System.currentTimeMillis()
def text = new File('/Users/pditommaso/Downloads/hs_ref_GRCh38_chr1.fa')
text.eachLine { String it ->
    data.put( bytes(count), KryoHelper.serialize(it) )
    count++
}

def end = System.currentTimeMillis()-now;
println "End - time: ${end/1000} secs"
data.close();
index.close()