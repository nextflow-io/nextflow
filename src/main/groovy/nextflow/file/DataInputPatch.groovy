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

package nextflow.file

import groovy.transform.CompileStatic
import org.mapdb.DataInput2

/**
 * Patch some bugs in MapDB {@link DataInput2} class
 * 1) DataInput2#read() do not return -1 when there's no more available data
 * 2) DataInput's methods do not throw {@link EOFException} when there's no more available data
 * 3) DataInput2#readUnsignedByte() move the read cursor forward even if there's no more available data
 *
 * See also https://github.com/jankotek/MapDB/pull/383
 *
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class DataInputPatch extends InputStream implements DataInput  {

    @Delegate
    DataInput2 target

    DataInputPatch(DataInput2 input) {
        this.target = input
    }

    protected void checkLimit(int i) {
        if ((i < 0) || (i >= target.buf.limit()))
            throw new EOFException();
    }

    protected int getPos() { target.pos }

    protected void setPos(int p) { target.pos = p }

    @Override
    public boolean readBoolean() throws IOException {
        checkLimit(pos)
        final result = target.buf.get(pos) == 1
        pos += 1
        return result
    }

    @Override
    public byte readByte() throws IOException {
        checkLimit(pos)
        final result = target.buf.get(pos)
        pos += 1
        return result
    }

    @Override
    public int readUnsignedByte() throws IOException {
        checkLimit(pos)
        final result = target.buf.get(pos) & 0xff;
        pos += 1
        return result
    }

    @Override
    public short readShort() throws IOException {
        checkLimit(pos)
        final result = target.buf.getShort(pos);
        pos += 2
        return result
    }

    @Override
    public int readUnsignedShort() throws IOException {
        checkLimit(pos)
        def result = (( (target.buf.get(pos) & 0xff) << 8) | ( (target.buf.get(pos+1) & 0xff)))
        pos += 2
        return result
    }

    @Override
    public int readInt() throws IOException {
        checkLimit(pos)
        final int ret = target.buf.getInt(pos);
        pos += 4
        return ret
    }

    @Override
    public long readLong() throws IOException {
        checkLimit(pos)
        final long ret = target.buf.getLong(pos);
        pos+=8;
        return ret;
    }

    @Override
    public float readFloat() throws IOException {
        checkLimit(pos)
        final float ret = target.buf.getFloat(pos);
        pos+=4;
        return ret;
    }

    @Override
    public double readDouble() throws IOException {
        checkLimit(pos)
        final double ret = target.buf.getDouble(pos);
        pos+=8;
        return ret;
    }


    @Override
    public int read() throws IOException {
        try {
            return readUnsignedByte()
        }
        catch( EOFException e ) {
            return -1
        }
    }

}
