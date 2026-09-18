/*
 * Copyright 2013-2026, Seqera Labs
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

package nextflow.cli

import groovy.transform.CompileStatic

/**
 * Read the user input from the console, falling back to the standard input
 * when no console is available e.g. the input is piped or redirected.
 *
 * The same instance must be shared by all the prompts of an interactive
 * session, since the underlying reader is stateful.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class ConsoleInput {

    private BufferedReader reader

    /**
     * @return The next line of user input, or {@code null} when no more input is available
     */
    String readLine() {
        final console = System.console()
        if( console != null )
            return console.readLine()
        // the reader is created once and reused for all the prompts, because it reads ahead
        // into an internal buffer -- a new reader for each line would drop the input that has
        // already been buffered, and therefore only the first line could be read
        if( reader == null )
            reader = new BufferedReader(new InputStreamReader(System.in))
        return reader.readLine()
    }

}
