/*
 * Copyright (c) 2013-2014, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2014, Paolo Di Tommaso and the respective authors.
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

package nextflow.trace
import java.text.DateFormat
import java.text.SimpleDateFormat

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.Memoized
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.util.Duration
import nextflow.util.MemoryUnit
/**
  * This object represent holds the information of a single process run,
  * its content is saved to a trace file line
  *
  * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
  */
@Slf4j
@CompileStatic
@EqualsAndHashCode(includes = 'store')
class TraceRecord {

    final private static String DEFAULT_DATE_FORMAT = "yyyy-MM-dd HH:mm:ss.SSS"

    final static NA = '-'

    final static List<String> NON_PRIMITIVE_TYPES = ['date','time','perc','mem']

    static final Map<String,String> FIELDS = [
            task_id:    'str',
            hash:       'str',
            native_id:  'str',
            name:       'str',
            status:     'str',
            exit_status:'num',
            submit:     'date',
            start:      'date',
            complete:   'date',
            wall_time:  'time',
            run_time:   'time',
            '%cpu':     'perc',     // -- ps field '%cpu'
            '%mem':     'perc',     // -- ps field '%mem'
            rss:        'mem',      // -- ps field 'rss'
            vmem:       'mem',      // -- ps field 'vsize'
            peak_rss:   'mem',      // -- /proc/$pid/status field 'VmHWM'  (Peak resident set size i.e. high water mark)
            peak_vmem:  'mem',      // -- /proc/$pid/status field 'VmPeak' (Peak virtual memory size)
            rchar:      'mem',      // -- /proc/$pid/io
            wchar:      'mem',      // -- /proc/$pid/io
            syscr:      'mem',      // -- /proc/$pid/io
            syscw:      'mem',      // -- /proc/$pid/io
            read_bytes: 'mem',      // -- /proc/$pid/io
            write_bytes:'mem'       // -- /proc/$pid/io
    ]

    static Map<String,Closure<String>> FORMATTER = [
            str: this.&fmtString,
            num: this.&fmtNumber,
            date: this.&fmtDate,
            time: this.&fmtTime,
            mem: this.&fmtMemory,
            perc: this.&fmtPercent
    ]

    @Memoized
    static private ThreadLocal<DateFormat> getLocalDateFormat(String fmt) {
        return new ThreadLocal<DateFormat>() {
            @Override
            protected DateFormat initialValue() {
                new SimpleDateFormat(fmt)
            }
        }
    }

    @PackageScope
    static DateFormat getDateFormat(String fmt = null) {
        if( !fmt )
            fmt = DEFAULT_DATE_FORMAT
        getLocalDateFormat(fmt).get()
    }


    /**
     * Convert the given value to a string
     *
     * @param fmt
     * @return The value as a string
     */
    @PackageScope
    static String fmtString(value, String fmt) {
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

        if( value instanceof Number )
            return value != Integer.MAX_VALUE ? value.toString() : NA

        if( value instanceof Duration )
            return value.toMillis().toString()

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
        getDateFormat(fmt).format(new Date(value as long))
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
    static def fmtPercent( def value, String fmt ) {
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
    Map<String,Object> store = [:]

    @Memoized
    def Set<String> keySet() {
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
    String get( String name, String converter ) {
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

    def getTaskId() { get('task_id') }

    /**
     * Render the specified list of fields to a single string value
     *
     * @param fields The list of fields to be rendered, each entry can optionally specify a
     *      a format string separating it from the name by a colon character e.g. name:format
     * @param separator A delimiter that separates fields entries in the final string
     * @return The final string containing all field values
     */
    String render( List<String> fields, List<String> formats, String separator ) {
        def result = new ArrayList(fields.size())
        for( int i=0; i<fields.size(); i++ ) {
            String name = fields[i]
            String format = i<formats?.size() ? formats[i] : null
            result << (get(name, format) ?: NA)
        }

        return result.join(separator)
    }


    String toString() {
        "${this.class.simpleName}${store}"
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

    TraceRecord parseTraceFile( String text ) {

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
                        this.put(name, values[i].toInteger() / 10F)
                    }
                    else if( i>3 ) {
                        def val = values[i].toLong()
                        // fields from index 4 to 7 (vmem,rss,peak_vmem, peak_rss) are provided in KB, so they are normalized to bytes
                        if( i<8 ) val *= 1024
                        this.put(name, val)
                    }

                }
            }

            // third line is supposed to be
            else if( count == 2 ) {
                try {
                    def elapsed = row.toString().trim().toLong()
                    this.put('run_time', elapsed)
                }
                catch( Exception e ) {
                    log.debug "Not a valid trace `run_time` value: '$count'"
                }
            }

        }

        return this
    }

}