/*
 * Copyright 2024, Ben Sherman
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

import java.nio.file.Path

import groovy.transform.CompileStatic

@CompileStatic
class CsvWriter {

    private /* boolean | List<String> */ header = false

    private String sep = ','

    CsvWriter(Map opts) {
        if( opts.header )
            this.header = opts.header

        if( opts.sep )
            this.sep = opts.sep.toString()
    }

    void apply(List records, Path path) {
        Collection columns
        if( header == true ) {
            final first = records.first()
            if( first !instanceof Map )
                throw new IllegalArgumentException('Records must be map objects when header=true')
            columns = ((Map)first).keySet()
        }
        else if( header instanceof List ) {
            columns = header
        }

        path.delete()

        if( columns )
            path << columns.collect(it -> '"' + it + '"').join(sep) << '\n'

        for( final record : records ) {
            Collection values
            if( record instanceof List )
                values = record
            else if( record instanceof Map )
                values = columns
                    ? record.subMap(columns).values()
                    : record.values()
            else
                throw new IllegalArgumentException('Records must be list or map objects')

            path << values.collect(it -> '"' + it + '"').join(sep) << '\n'
        }
    }

}
