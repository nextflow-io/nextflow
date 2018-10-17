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

package nextflow.io

import groovy.transform.CompileStatic

/**
 * Implements a  {@link InputStream} object backed on a {@link DataInput} instance
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class DataInputStreamAdapter extends InputStream {

    final private DataInput target

    DataInputStreamAdapter( DataInput input ) { this.target = input }

    @Override
    int read() throws IOException {
        try {
            return target.readUnsignedByte()
        }
        catch( EOFException | IndexOutOfBoundsException e ) {
            return -1
        }
    }

    long skip(long n) {
        target.skipBytes(n as int)
    }
}
