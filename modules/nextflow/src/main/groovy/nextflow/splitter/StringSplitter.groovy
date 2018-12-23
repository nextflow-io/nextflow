/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

package nextflow.splitter
import groovy.transform.CompileStatic
import groovy.transform.InheritConstructors
import groovy.util.logging.Slf4j
/**
 * Simple slitter chunking a string in sub-strings having the specified length
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
@InheritConstructors
class StringSplitter extends AbstractTextSplitter {

    protected boolean ignoreNewLine

    StringSplitter options(Map options) {
        super.options(options)
        ignoreNewLine = options.ignoreNewLine == true ?: false
        return this
    }

    /**
     * @return A map representing the valid options for the splitter. The map keys define the
     * accepted parameter names, the values the valid values for each of them.
     */
    @Override
    protected Map<String,Object> validOptions() {
        def result = super.validOptions()
        result.ignoreNewLine = Boolean
        return result
    }


    @Override
    protected fetchRecord(BufferedReader targetObject) {

        int ch
        while( true ) {
            ch = targetObject.read()
            if( ch == -1 )
                return null

            if( ignoreNewLine && ( ch == '\n' as char || ch == '\r' as char ))
                continue

            break
        }

        return ch as char
    }
}
