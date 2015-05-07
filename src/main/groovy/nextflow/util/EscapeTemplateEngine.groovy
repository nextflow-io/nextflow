/*
 * Copyright (c) 2013-2015, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2015, Paolo Di Tommaso and the respective authors.
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

/*
 * Copyright 2003-2013 the original author or authors.
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

package nextflow.util

import groovy.text.Template
import groovy.text.TemplateEngine
import org.codehaus.groovy.control.CompilationFailedException
import org.codehaus.groovy.runtime.InvokerHelper

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
public class EscapeTemplateEngine extends TemplateEngine {
    private boolean verbose;
    private static int counter = 1;

    private GroovyShell groovyShell;

    public EscapeTemplateEngine() {
        this(GroovyShell.class.getClassLoader());
    }

    public EscapeTemplateEngine(boolean verbose) {
        this(GroovyShell.class.getClassLoader());
        setVerbose(verbose);
    }

    public EscapeTemplateEngine(ClassLoader parentLoader) {
        this(new GroovyShell(parentLoader));
    }

    public EscapeTemplateEngine(GroovyShell groovyShell) {
        this.groovyShell = groovyShell;
    }

    public Template createTemplate(Reader reader) throws CompilationFailedException, IOException {
        SimpleTemplate template = new SimpleTemplate();
        String script = template.parse(reader);
        if (verbose) {
            System.out.println("\n-- script source --");
            System.out.print(script);
            System.out.println("\n-- script end --\n");
        }
        try {
            template.script = groovyShell.parse(script, "EscapeTemplateScript" + counter++ + ".groovy");
        } catch (Exception e) {
            throw new GroovyRuntimeException("Failed to parse template script (your template may contain an error or be trying to use expressions not currently supported): " + e.getMessage());
        }
        return template;
    }

    /**
     * @param verbose true if you want the engine to display the template source file for debugging purposes
     */
    public void setVerbose(boolean verbose) {
        this.verbose = verbose;
    }

    public boolean isVerbose() {
        return verbose;
    }

    private static class SimpleTemplate implements Template {

        protected Script script;

        public Writable make() {
            return make(null);
        }

        final static private DOLLAR = (int)('$' as char)
        final static private BACKSLASH = (int)('\\' as char)
        final static private PLACEHOLDER = (int)('&' as char)
        final static private CURLY_OPEN = (int)('{' as char)
        final static private CURLY_CLOSE = (int)('}' as char)
        final static private DOUBLE_QUOTE = (int)('"' as char)
        final static private NL = (int)('\n' as char)
        final static private CR = (int)('\r' as char)
        final static private PERIOD = (int)('.' as char)


        public Writable make(final Map map) {
            return new Writable() {
                /**
                 * Write the template document with the set binding applied to the writer.
                 *
                 * @see groovy.lang.Writable#writeTo(java.io.Writer)
                 */
                public Writer writeTo(Writer writer) {
                    Binding binding;
                    if (map == null)
                        binding = new Binding();
                    else
                        binding = new Binding(map);
                    Script scriptObject = InvokerHelper.createScript(script.getClass(), binding);
                    PrintWriter pw = new PrintWriter(writer);
                    scriptObject.setProperty("out", pw);
                    scriptObject.run();
                    pw.flush();
                    return writer;
                }

                /**
                 * Convert the template and binding into a result String.
                 *
                 * @see java.lang.Object#toString()
                 */
                public String toString() {
                    StringWriter sw = new StringWriter();
                    writeTo(sw);
                    return sw.toString();
                }
            };
        }

        /**
         * Parse the text document looking for <% or <%= and then call out to the appropriate handler, otherwise copy the text directly
         * into the script while escaping quotes.
         *
         * @param reader a reader for the template text
         * @return the parsed text
         * @throws IOException if something goes wrong
         */
        protected String parse(Reader reader) throws IOException {
            if (!reader.markSupported()) {
                reader = new BufferedReader(reader);
            }
            StringWriter sw = new StringWriter();
            startScript(sw);
            int c;
            while ((c = reader.read()) != -1) {

                if (c == DOLLAR) {
                    // escape dollar characters
                    sw.write(BACKSLASH)
                    sw.write(DOLLAR)
                    continue
                }

                if (c == BACKSLASH) {
                    // escape backslash itself
                    sw.write(BACKSLASH)
                    sw.write(BACKSLASH)
                    continue
                }

                if (c == PLACEHOLDER) {
                    reader.mark(1);
                    c = reader.read();      // read the next character
                    if ( c == CURLY_OPEN ){
                        reader.mark(1)
                        sw.write(DOLLAR)    // <-- replace the placeholder with a $ char
                        sw.write(CURLY_OPEN)
                        processGSstring(reader, sw);
                    }
                    else if (Character.isJavaIdentifierStart(c)) {
                        sw.write(DOLLAR)    // <-- replace the placeholder with a $ char
                        reader.reset()
                        processIdentifier(reader, sw)
                    }
                    else {
                        // just write this char and continue
                        sw.write(PLACEHOLDER)
                        sw.write(c)
                        // escape back slash doubling it
                        if (c == BACKSLASH) sw.write(c)
                    }

                    continue; // at least '$' is consumed ... read next chars.
                }

                if (c == DOUBLE_QUOTE) {
                    sw.write('\\');
                }

                /*
                 * Handle raw new line characters.
                 */
                if (c == NL || c == CR) {
                    if (c == CR) { // on Windows, "\r\n" is a new line.
                        reader.mark(1);
                        c = reader.read();
                        if (c != NL) {
                            reader.reset();
                        }
                    }
                    sw.write("\n");
                    continue;
                }
                sw.write(c);
            }
            endScript(sw);
            return sw.toString();
        }

        private void startScript(StringWriter sw) {
            sw.write("out.print(\"\"\"");
        }

        private void endScript(StringWriter sw) {
            sw.write("\"\"\");\n");
            sw.write("\n/* Generated by SimpleTemplateEngine */");
        }

        private void processGSstring(Reader reader, StringWriter sw) throws IOException {
            int c;
            def name = new StringBuilder()

            while ((c = reader.read()) != -1) {
                if (c != NL && c != CR) {
                    sw.write(c);
                }
                if (c == CURLY_CLOSE) {
                    break;
                }
                name.append((char)c)
            }

        }

        private void processIdentifier(Reader reader, StringWriter sw) {
            int c;
            int pos=0;
            def name = new StringBuilder()

            while ((c = reader.read()) != -1) {
                sw.write(c);

                if( (pos==0 && Character.isJavaIdentifierStart(c)) || Character.isJavaIdentifierPart(c) ) {
                    pos++
                    name.append((char)c)
                    continue
                }

                if( pos && c == PERIOD ) {
                    reader.mark(1)
                    char next = reader.read()
                    if( Character.isJavaIdentifierStart(next)) {
                        pos = 0;
                        sw.write(next)
                        name.append('.')
                        name.append(next)
                        continue
                    }
                    else {
                        reader.reset()
                    }
                }
                break
            }

        }


    }
}