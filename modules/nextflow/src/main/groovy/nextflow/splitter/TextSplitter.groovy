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

package nextflow.splitter

import groovy.transform.CompileStatic
import groovy.transform.InheritConstructors
import groovy.util.logging.Slf4j
/**
 * Split a text file by one or more lines at times
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
@InheritConstructors
class TextSplitter extends AbstractTextSplitter {

    private boolean keepHeader

    @Override
    protected Map<String,Object> validOptions() {
        def result = super.validOptions()
        result.keepHeader = [Boolean]
        return result
    }

    @Override
    TextSplitter options(Map opts) {
        super.options(opts)
        this.keepHeader = opts.keepHeader == true
        return this
    }

    /**
     * A record is a text line
     *
     * @param reader The buffered reader
     * @return A line string or {@code null} when the end of the file is reached
     */
    @Override
    protected fetchRecord(BufferedReader reader) {
        def line = reader.readLine()
        if( line != null ) line+='\n'
        return line
    }

    protected void parseHeader(Reader reader) {
        if( !keepHeader )
            return
        
        def line = reader.readLine()
        if( line==null )
            return

        def collector = getCollector()
        if( collector instanceof HeaderCollector ) {
            collector.setHeader(line+'\n')
        }
    }

    @Override
    protected process( Reader targetObject ) {
        final reader = wrapReader(targetObject)
        parseHeader(reader)
        super.process(reader)
    }
}
