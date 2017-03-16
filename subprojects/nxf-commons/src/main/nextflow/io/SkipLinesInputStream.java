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

    public long skip(long n) throws IOException {
        throw new UnsupportedOperationException("Skip operation is not supported by " + this.getClass().getName());
    }

    @Override
    public int read() throws IOException {

        while( count < skip ) {
            int ch  = read0();
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

        int result = target.read();
        if( result == LF ) {
            count++;
        }
        if( result == CR ) {
            count++;
            int next = target.read();
            if( next != LF )
                buffer = next;
        }
        return result;
    }
}
