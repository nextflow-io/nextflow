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
        path.delete()

        final columns = columnHeaders(header, records)
        if( columns )
            path << columns.collect(column -> "\"${column}\"").join(sep) << '\n'

        if( records.isEmpty() )
            path << ''

        for( final record : records ) {
            final values = rowValues(record, columns)
            path << values.collect(v -> "\"${toCsvString(v)}\"").join(sep) << '\n'
        }
    }

    private static Collection columnHeaders(Object header, List records) {
        if( header instanceof List ) {
            return header
        }

        if( header == true && !records.isEmpty() ) {
            final first = records.first()
            if( first instanceof Map )
                return first.keySet()
            else
                throw new IllegalArgumentException('Records must be map objects when header=true')
        }

        return null
    }

    private static Collection rowValues(Object record, Collection columns) {
        if( record instanceof List ) {
            return record
        }

        if( record instanceof Map ) {
            return columns
                ? record.subMap(columns).values()
                : record.values()
        }

        if( isSerializable(record) ) {
            return [ record ]
        }

        throw new IllegalArgumentException("Record of type `${record.class.name}` can not be serialized to CSV")
    }

    private static boolean isSerializable(value) {
        return value == null
            || value instanceof Boolean
            || value instanceof CharSequence
            || value instanceof Number
            || value instanceof Path
    }

    private static String toCsvString(value) {
        if( value == null )
            return ""

        if( value instanceof Path )
            return value.toUriString()

        return value.toString()
    }

}
