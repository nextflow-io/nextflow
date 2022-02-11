/*
 * Copyright 2020-2022, Seqera Labs
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/*
 * Copyright 2012-2014 the original author or authors.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package test;

import org.junit.rules.TestRule;
import org.junit.runner.Description;
import org.junit.runners.model.Statement;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;

/**
 * JUnit {@code @Rule} to capture output from System.out and System.err.
 *
 * @author Phillip Webb
 */
public class OutputCapture implements TestRule {

    private CaptureOutputStream captureOut;

    private CaptureOutputStream captureErr;

    private ByteArrayOutputStream copy;

    @Override
    public Statement apply(final Statement base, Description description) {
        return new Statement() {
            @Override
            public void evaluate() throws Throwable {
                captureOutput();
                try {
                    base.evaluate();
                }
                finally {
                    releaseOutput();
                }
            }
        };
    }

    protected void captureOutput() {
        AnsiOutputControl.get().disableAnsiOutput();
        this.copy = new ByteArrayOutputStream();
        this.captureOut = new CaptureOutputStream(System.out, this.copy);
        this.captureErr = new CaptureOutputStream(System.err, this.copy);
        System.setOut(new PrintStream(this.captureOut));
        System.setErr(new PrintStream(this.captureErr));
    }

    protected void releaseOutput() {
        AnsiOutputControl.get().enabledAnsiOutput();
        System.setOut(this.captureOut.getOriginal());
        System.setErr(this.captureErr.getOriginal());
        this.copy = null;
    }

    public void flush() {
        try {
            this.captureOut.flush();
            this.captureErr.flush();
        }
        catch (IOException ex) {
            // ignore
        }
    }

    @Override
    public String toString() {
        flush();
        return this.copy.toString();
    }

    private static class CaptureOutputStream extends OutputStream {

        private final PrintStream original;

        private final OutputStream copy;

        public CaptureOutputStream(PrintStream original, OutputStream copy) {
            this.original = original;
            this.copy = copy;
        }

        @Override
        public void write(int b) throws IOException {
            this.copy.write(b);
            this.original.write(b);
            this.original.flush();
        }

        @Override
        public void write(byte[] b) throws IOException {
            write(b, 0, b.length);
        }

        @Override
        public void write(byte[] b, int off, int len) throws IOException {
            this.copy.write(b, off, len);
            this.original.write(b, off, len);
        }

        public PrintStream getOriginal() {
            return this.original;
        }

        @Override
        public void flush() throws IOException {
            this.copy.flush();
            this.original.flush();
        }

    }

    /**
     * Allow AnsiOutput to not be on the test classpath.
     */
    private static class AnsiOutputControl {

        public void disableAnsiOutput() {
        }

        public void enabledAnsiOutput() {
        }

        public static AnsiOutputControl get() {
            return new AnsiOutputControl();
        }

    }



}
