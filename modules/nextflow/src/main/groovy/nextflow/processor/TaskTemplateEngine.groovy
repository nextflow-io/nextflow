/*
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

package nextflow.processor
import ch.grengine.Grengine
import groovy.text.Template
import groovy.text.TemplateEngine
import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import org.codehaus.groovy.control.CompilationFailedException
/**
 * A template engine that uses {@link Grengine} to parse template scripts
 * It also that ignore dollar variable and uses a custom character for variable interpolation
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class TaskTemplateEngine extends TemplateEngine {

    final static private int DOLLAR = (int)('$' as char)
    final static private int BACKSLASH = (int)('\\' as char)
    final static private int CURLY_OPEN = (int)('{' as char)
    final static private int CURLY_CLOSE = (int)('}' as char)
    final static private int DOUBLE_QUOTE = (int)('"' as char)
    final static private int NL = (int)('\n' as char)
    final static private int CR = (int)('\r' as char)
    final static private int PERIOD = (int)('.' as char)

    private Grengine grengine;

    private char placeholder = '$' as char

    private boolean enableShortNotation

    private Set<String> variableNames

    private String result

    TaskTemplateEngine() {
        grengine = new Grengine()
    }

    TaskTemplateEngine( Grengine engine ) {
        this.grengine = engine
    }

    String getResult() { result }

    Set<String> getVariableNames() { variableNames }

    TaskTemplateEngine eval(String text, Map binding = null) {
        def template = createTemplate(text)
        def result = binding ? template.make(binding) : template.make()
        this.result = result?.toString()
        return this
    }

    String render( String text, Map binding = null )  {
        def template = createTemplate(text)
        def result = binding ? template.make(binding) : template.make()
        result?.toString()
    }

    Template createTemplate(Reader reader) throws CompilationFailedException, IOException {
        ParsableTemplate template = placeholder == DOLLAR ? new SimpleTemplate() : new EscapeTemplate()
        template.create(reader, grengine)
        variableNames = template.getVariablesNames()
        return template;
    }

    TaskTemplateEngine setEnableShortNotation(boolean value) {
        enableShortNotation = value
        return this
    }

    TaskTemplateEngine setPlaceholder(char ch) {
        placeholder = ch
        return this
    }

    /**
     * Common template methods
     */
    @Slf4j
    private static abstract class ParsableTemplate implements Template {

        protected Script script

        protected Grengine grengine

        Set<String> getVariablesNames() {
            def result = script.getBinding().getVariables().get('__$$_template_vars') as Set<String>
            return result
        }

        abstract String parse(Reader reader) throws IOException

        void create(Reader reader, Grengine grengine) {

            // -- set the grengine instance
            this.grengine = grengine

            // -- parse the template reader and create script string
            String script = parse(reader);
            log.trace "\n-- script source --\n${script}\n-- script end --\n"

            // -- finally create the Script object
            try {
                this.script = grengine.create(script);
            }
            catch (Exception e) {
                throw new GroovyRuntimeException("Failed to parse template script (your template may contain an error or be trying to use expressions not currently supported): " + e.getCause() ?: e.toString());
            }

        }

        protected void startScript(StringWriter sw) {
            sw.write('__$$_out.print("""');
        }

        protected void endScript(StringWriter sw) {
            sw.write('""");\n');
            sw.write('\n');
        }


        Writable make() {
            return make(null);
        }

        Writable make(final Map context) {
            return new Writable() {

                /**
                 * Write the template document with the set binding applied to the writer.
                 *
                 * @see groovy.lang.Writable#writeTo(java.io.Writer)
                 */
                Writer writeTo(Writer writer) {

                    final binding = new TemplateBinding(context,writer)
                    try {
                        grengine.run(script, binding)
                    }
                    finally {
                        binding.flush();
                    }
                    return writer;
                }

                /**
                 * Convert the template and binding into a result String.
                 *
                 * @see java.lang.Object#toString()
                 */
                String toString() {
                    StringWriter sw = new StringWriter();
                    writeTo(sw);
                    return sw.toString();
                }
            }
        }
    }

    /**
     * Extends the {@link Binding} class adding an implicit writer object,
     * named {@code __$$_out} that is used by the template engine to render the target script,
     * in order to avoid to change the original context map
     */
    @PackageScope
    static class TemplateBinding extends Binding implements Closeable {

        private PrintWriter writer

        TemplateBinding(Map context, Writer writer) {
            super(context != null ? context : [:])
            this.writer = writer instanceof PrintWriter ? (PrintWriter)writer : new PrintWriter(writer)
        }

        @Override
        Object getVariable(String name) {
            if( name == '__$$_out' )
                return writer
            else
                return super.getVariable(name)
        }

        final void flush() {
            writer.flush()
        }

        void close() throws IOException {
            writer.close()
        }

    }

    /**
     * Template class escaping standard $ prefixed variables and using a custom character
     * as variable placeholder
     */
    private class EscapeTemplate extends ParsableTemplate {

        /**
         * Parse the text document looking for <% or <%= and then call out to the appropriate handler, otherwise copy the text directly
         * into the script while escaping quotes.
         *
         * @param reader a reader for the template text
         * @return the parsed text
         * @throws IOException if something goes wrong
         */
        String parse(Reader reader) throws IOException {
            if (!reader.markSupported()) {
                reader = new BufferedReader(reader);
            }

            final PLACEHOLDER = (int)TaskTemplateEngine.this.placeholder

            StringWriter sw = new StringWriter();
            startScript(sw);
            int c;
            while ((c = reader.read()) != -1) {

                if( c == DOLLAR ) {
                    // escape dollar characters
                    sw.write(BACKSLASH)
                    sw.write(DOLLAR)
                    continue
                }

                if( c == BACKSLASH ) {
                    // escape backslash itself
                    sw.write(BACKSLASH)
                    sw.write(BACKSLASH)
                    continue
                }

                if( c == PLACEHOLDER ) {
                    reader.mark(1);
                    c = reader.read();      // read the next character
                    if( c == CURLY_OPEN ){
                        reader.mark(1)
                        sw.write(DOLLAR)    // <-- replace the placeholder with a $ char
                        sw.write(CURLY_OPEN)
                        processGString(reader, sw);
                    }
                    else if( enableShortNotation && Character.isJavaIdentifierStart(c) ) {
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

                if( c == DOUBLE_QUOTE ) {
                    sw.write('\\');
                }

                /*
                 * Handle raw new line characters.
                 */
                if( c == NL || c == CR ) {
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

        private void processGString(Reader reader, StringWriter sw) throws IOException {
            int c
            while ((c = reader.read()) != -1) {
                if( c != NL && c != CR ) {
                    sw.write(c);
                }
                if(c == CURLY_CLOSE ) {
                    break;
                }
            }
        }

        private void processIdentifier(Reader reader, StringWriter sw) {
            int c
            int pos=0

            while ((c = reader.read()) != -1) {
                sw.write(c);

                if( (pos==0 && Character.isJavaIdentifierStart(c)) || Character.isJavaIdentifierPart(c) ) {
                    pos++
                    continue
                }

                if( pos && c == PERIOD ) {
                    reader.mark(1)
                    char next = reader.read()
                    if( Character.isJavaIdentifierStart(next)) {
                        pos = 0;
                        sw.write(next)
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

    /**
     * Default template that interpolates variables
     */
    private static class SimpleTemplate extends ParsableTemplate {

        /**
         * Parse the text document looking for <% or <%= and then call out to the appropriate handler, otherwise copy the text directly
         * into the script while escaping quotes.
         *
         * @param reader a reader for the template text
         * @return the parsed text
         * @throws IOException if something goes wrong
         */
        String parse(Reader reader) throws IOException {
            if (!reader.markSupported()) {
                reader = new BufferedReader(reader);
            }
            StringWriter sw = new StringWriter();
            startScript(sw);
            int c;
            while ((c = reader.read()) != -1) {

                if (c == DOLLAR) {
                    reader.mark(1);
                    c = reader.read();
                    if (c != CURLY_OPEN) {
                        sw.write(DOLLAR);
                        reader.reset();
                    } else {
                        reader.mark(1);
                        sw.write(DOLLAR);
                        sw.write(CURLY_OPEN);

                        processGString(reader, sw);
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
                        if (c != '\n') {
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


        private void processGString(Reader reader, StringWriter sw) throws IOException {
            int c;
            while ((c = reader.read()) != -1) {
                if (c != NL && c != CR) {
                    sw.write(c);
                }
                if (c == CURLY_CLOSE) {
                    break;
                }
            }
        }
    }

}