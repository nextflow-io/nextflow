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

package nextflow.trace
import java.nio.file.Path

import groovy.json.StringEscapeUtils
import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.extension.Bolts
import nextflow.processor.TaskId
import nextflow.util.Duration
import nextflow.util.KryoHelper
import nextflow.util.MemoryUnit
/**
  * This object represent holds the information of a single process run,
  * its content is saved to a trace file line
  *
  * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
  */
@Slf4j
@CompileStatic
class TraceRecord implements Serializable {

    TraceRecord() {
        this.store = [:]
    }

    private TraceRecord(Map store) {
        this.store = store
    }

    final private static String DEFAULT_DATE_FORMAT = "yyyy-MM-dd HH:mm:ss.SSS"

    final public static String NA = '-'

    final public static List<String> NON_PRIMITIVE_TYPES = ['date','time','perc','mem']

    final public static Map<String,String> FIELDS = [
            task_id:    'str',
            hash:       'str',
            native_id:  'str',
            process:    'str',
            module:     'str',
            container:  'str',
            tag:        'str',
            name:       'str',
            status:     'str',
            exit:       'str',
            submit:     'date',
            start:      'date',
            complete:   'date',
            duration:   'time',
            realtime:   'time',
            '%cpu':     'perc',     // -- ps field '%cpu'
            '%mem':     'perc',     // -- ps field '%mem'
            rss:        'mem',      // -- ps field 'rss'
            vmem:       'mem',      // -- ps field 'vsize'
            peak_rss:   'mem',      // -- /proc/$pid/status field 'VmHWM'  (Peak resident set size i.e. high water mark)
            peak_vmem:  'mem',      // -- /proc/$pid/status field 'VmPeak' (Peak virtual memory size)
            rchar:      'mem',      // -- /proc/$pid/io
            wchar:      'mem',      // -- /proc/$pid/io
            syscr:      'num',      // -- /proc/$pid/io
            syscw:      'num',      // -- /proc/$pid/io
            read_bytes: 'mem',      // -- /proc/$pid/io
            write_bytes:'mem',      // -- /proc/$pid/io
            attempt:    'num',
            workdir:    'str',
            script:     'str',
            scratch:    'str',
            queue:      'str',
            cpus:       'num',
            memory:     'mem',
            disk:       'mem',
            time:       'time',
            env:        'str',
            error_action:'str'
    ]

    static public Map<String,Closure<String>> FORMATTER = [
            str: this.&fmtString,
            num: this.&fmtNumber,
            date: this.&fmtDate,
            time: this.&fmtTime,
            mem: this.&fmtMemory,
            perc: this.&fmtPercent
    ]

    @PackageScope
    static TimeZone TIMEZONE = null


    /**
     * Convert the given value to a string
     *
     * @param fmt
     * @return The value as a string
     */
    @PackageScope
    static String fmtString(value, String fmt) {
        if( value instanceof Number )
            return value != Integer.MAX_VALUE ? value.toString() : NA

        value ? value.toString() : NA
    }


    /**
     * Convert the given value to a number.
     * <p>Objects of type {@link Duration} (time) are  converted to milliseconds
     * <p>Objects of type {@link MemoryUnit) are converted to number of bytes
     * <p>Objects of type {@link Date) are converted to milliseconds since Unix epoch
     *
     * @param fmt
     * @return The value as a string
     */
    @PackageScope
    static String fmtNumber(def value, String fmt) {
        if( value == null )
            return NA

       if( value instanceof Duration )
            return String.valueOf(value.toMillis())

        if( value instanceof MemoryUnit )
            return String.valueOf(value.toBytes())

        if( value instanceof Date )
            return String.valueOf(value.getTime())

        return value.toString()
    }

    /**
     * Converts the value to a date string
     *
     * @param value The value is supposed to be the number of milliseconds since Unix epoch
     * @param fmt
     * @return The formatted date string
     */
    @PackageScope
    static String fmtDate(def value, String fmt) {
        if( !value ) return NA
        if( !fmt )
            fmt = DEFAULT_DATE_FORMAT
        Bolts.format(new Date(value as long), fmt, TIMEZONE)
    }

    /**
     * Coverts the value to a duration string.
     *
     * See {@link Duration}
     * @param value
     * @param fmt
     * @return
     */
    @PackageScope
    static String fmtTime(def value, String fmt) {
        if( value == null ) return NA
        new Duration(value as long).toString()
    }

    /**
     * Converts the value to a memory unit string
     * <p>
     * See {@link MemoryUnit}
     *
     * @param value
     * @param fmt
     * @return
     */
    @PackageScope
    static String fmtMemory( def value, String fmt) {
        if( value == null ) return NA

        String str = value.toString()
        if( str.isLong() ) {
            str = new MemoryUnit(str.toLong()).toString()
            str = str.replaceAll(/,/,'')
        }

        return str
    }

    /**
     * Converts the value to a percent value string
     * @param value
     * @param fmt
     * @return
     */
    @PackageScope
    static String fmtPercent( def value, String fmt ) {
        if( value == null ) return NA
        try {
            if( value instanceof Number )
                return String.format('%.1f%%', value.toFloat())
            else {
                def x = value.toString().toFloat()
                return String.format('%.1f%%', x)
            }

        }
        catch( Exception e ) {
            log.debug "Not a valid percentual value: '$value'"
            return NA
        }
    }


    @PackageScope
    Map<String,Object> store

    @Memoized
    Set<String> keySet() {
        FIELDS.keySet()
    }

    def propertyMissing(String name, value) {
        put(name,value)
    }

    def propertyMissing(String name) {
        get(name)
    }

    def containsKey( String name ) {
        assert name
        store.containsKey(name)
    }

    def get( String name ) {
        assert keySet().contains(name), "Not a valid TraceRecord field: '$name'"
        store.get(name)
    }

    def void put( String name, def value ) {
        assert keySet().contains(name), "Not a valid TraceRecord field: '$name'"

        // vmpeak: Peak virtual memory size
        // this is a synonym of 'max_vmem' field
        if( name == 'vmpeak' ) {
            store.put('max_vmem', value)
        }

        // Peak resident set size ("high water mark")
        // This is a synonym of 'max_rss' field
        else if( name == 'vmhwm' ) {
            store.put('max_rss', value)
        }

        else {
            store.put(name, value)
        }
    }

    def void putAll( Map<String,Object> values ) {
        if( !values )
            return

        for( String key : values.keySet() ) {
            put(key, values.get(key))
        }
    }

    /**
     * Get a trace field value and apply a conversion rule to it
     *
     * @param name The field name e.g. task_id, status, etc.
     * @param converter A converter string
     * @return A string value formatted according the specified converter
     */
    String getFmtStr( String name, String converter = null ) {
        assert name
        def val = store.get(name)

        String sType=null
        String sFormat=null
        if( converter ) {
            int p = converter.indexOf(':')
            if( p == -1 ) {
                sType = converter
            }
            else {
                sType = converter.substring(0,p)
                sFormat = converter.substring(p+1)
            }
        }

        def type = sType ?: FIELDS.get(name)
        if( !type )
            throw new IllegalArgumentException("Not a valid trace field name: '$name'")


        def formatter = FORMATTER.get(type)
        if( !formatter )
            throw new IllegalArgumentException("Not a valid trace formatter for field: '$name' with type: '$type'")

        try {
            return formatter.call(val,sFormat)
        }
        catch( Throwable e ) {
            log.debug "Not a valid trace value -- field: '$name'; value: '$val'; format: '$sFormat'"
            return null
        }
    }

    TaskId getTaskId() { (TaskId)get('task_id') }

    String getWorkDir() { get('workdir') }

    /**
     * Render the specified list of fields to a single string value
     *
     * @param fields The list of fields to be rendered, each entry can optionally specify a
     *      a format string separating it from the name by a colon character e.g. name:format
     * @param separator A delimiter that separates fields entries in the final string
     * @return The final string containing all field values
     */
    String renderText( List<String> fields, List<String> formats, String separator ) {
        def result = new ArrayList(fields.size())
        for( int i=0; i<fields.size(); i++ ) {
            String name = fields[i]
            String format = i<formats?.size() ? formats[i] : null
            result << (getFmtStr(name, format) ?: NA)
        }

        return result.join(separator)
    }

    CharSequence renderJson(StringBuilder result, List<String> fields, List<String> formats) {
        final QUOTE = '"'
        if( result == null ) result = new StringBuilder()
        result << "{"
        for( int i=0; i<fields.size(); i++ ) {
            if( i ) result << ','
            String name = fields[i]
            String format = i<formats?.size() ? formats[i] : null
            String value = StringEscapeUtils.escapeJavaScript(getFmtStr(name, format) ?: NA)
            result << QUOTE << name << QUOTE << ":" << QUOTE << value << QUOTE
        }
        result << "}"
        return result
    }

    CharSequence renderJson(StringBuilder result = new StringBuilder()) {
        def fields = []
        def formats = []
        FIELDS.each { name, type -> fields << name; formats << type }
        renderJson(result, fields, formats)
    }

    String toString() {
        "${this.class.simpleName} ${store}"
    }


    /**
     * Parse the trace file
     * <p>
     * Trace example:
     * <pre>
     * pid state %cpu %mem vmem rss max_vmem max_rss rchar wchar syscr syscw read_bytes write_bytes
     *  1 0 0 0 11084 1220 11084 1220 4790 12 11 1 0 0 0
     * </pre>
     *
     */

    TraceRecord parseTraceFile( Path file ) {

        final text = file.text
        String[] header = null
        final lines = text.readLines()
        for( int count=0; count<lines.size(); count++ ) {
            String row = lines[count]

            /*
             * 1st line -- parse the header
             */
            if( count == 0 ) {
                header = row.trim().split(/\s+/)
            }

            /*
             * 2nd line -- parse values produced by 'ps'
             */
            else if( count == 1 ) {
                String[] values = row.trim().split(/\s+/)
                for( int i=0; i<values.length; i++ ) {

                    final name = header[i]
                    if( i==2 || i==3 ) {
                        // fields '%cpu' and '%mem' are expressed as percent value
                        this.put(name, parseInt(values[i], file, i) / 10F)
                    }
                    else if( i>3 ) {
                        def val = parseLong(values[i], file, i)
                        // fields from index 4 to 7 (vmem,rss,peak_vmem, peak_rss) are provided in KB, so they are normalized to bytes
                        if( i<8 ) val *= 1024
                        this.put(name, val)
                    }

                }
            }

            // third line is the cpu realtime
            else if( count == 2 ) {
                try {
                    def elapsed = row.toString().trim().toLong()
                    this.put('realtime', elapsed)
                }
                catch( Exception e ) {
                    log.debug "Not a valid trace `realtime` value: '$count'"
                }
            }

        }

        return this
    }

    private long parseInt( String str, Path file, int index )  {
        try {
            str.toInteger()
        }
        catch( NumberFormatException e ) {
            log.debug "[WARN] Not a valid integer number `$str` -- offending column: $index in file `$file`"
            return 0
        }
    }

    private long parseLong( String str, Path file, int index )  {
        try {
            str.toLong()
        }
        catch( NumberFormatException e ) {
            log.debug "[WARN] Not a valid long number `$str` -- offending column: $index in file `$file`"
            return 0
        }
    }

    @Override
    boolean equals( Object that ) {
        if ( this.is(that) ) return true
        if ( !(that instanceof TraceRecord) ) return false
        this.store.equals(((TraceRecord)that).store)
    }

    @Override
    int hashCode() {
        store.hashCode()
    }

    byte[] serialize() {
        KryoHelper.serialize(store)
    }

    static TraceRecord deserialize(byte[] buffer) {
        Map map = (Map)KryoHelper.deserialize(buffer)
        new TraceRecord(map)
    }

    TraceRecord setCached(boolean value) {
        if( value ) {
            store.status = 'CACHED'
        }
        return this
    }

    boolean isCached() {
        store.status == 'CACHED'
    }

}
