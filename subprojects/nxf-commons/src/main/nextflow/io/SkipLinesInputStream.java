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
package nextflow.io;

import java.io.IOException;
import java.io.InputStream;

/**
 * An input stream that skips the first `n` lines in a text stream
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
public class SkipLinesInputStream extends InputStream {

    private static final int LF = '\n';

    private static final int CR = '\r';

    private InputStream target;

    private int skip;

    private int buffer=-1;

    private int count;

    private StringBuilder header;

    /**
     * Creates a <code>FilterInputStream</code>
     * by assigning the  argument <code>in</code>
     * to the field <code>this.in</code> so as
     * to remember it for later use.
     *
     * @param inputStream the underlying input stream, or <code>null</code> if
     *          this instance is to be created without an underlying stream.
     */
    public SkipLinesInputStream(InputStream inputStream, int skip) {
        this.target = inputStream;
        this.skip = skip;
    }

    @Override
    public void close() throws IOException {
        target.close();
    }

    public long skip(long n) throws IOException {
        throw new UnsupportedOperationException("Skip operation is not supported by " + this.getClass().getName());
    }

    @Override
    public int read() throws IOException {

        while( count < skip ) {
            int ch = read0();
            if( ch == -1 )
                return -1;
        }

        return read0();
    }

    /**
     * A line is considered to be terminated by any one
     * of a line feed ('\n'), a carriage return ('\r'), or a carriage return
     * followed immediately by a linefeed.
     *
     */
    private int read0() throws IOException {
        if( buffer != -1 ) {
            int result = buffer;
            buffer=-1;
            return result;
        }

        int ch = target.read();
        if( header!=null && count<skip && ch!=-1 )
            header.append((char)ch);

        if( ch == LF ) {
            count++;
        }
        else if( ch == CR ) {
            count++;
            int next = target.read();
            if( next == LF ) {
                if( header!=null )
                    header.append((char)next);
            }
            else {
                buffer = next;
                if( header!=null && count<skip && next!=-1 )
                    header.append((char)next);
            }

        }
        return ch;
    }

    public String consumeHeader() throws IOException {
        header = new StringBuilder();
        while( count<skip && read0()!=-1 ) ;
        return header.toString();
    }

    public String getHeader() {
        return header != null ? header.toString() : null;
    }
}
