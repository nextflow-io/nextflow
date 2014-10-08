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
class TraceRecord {

    final private static String DEFAULT_DATE_FORMAT = "yyyy-MM-dd HH:mm:ss.SSS"

    final static NA = '-'

    static final Map<String,String> FIELDS = [
            task_id:    'str',
            hash:       'str',
            native_id:  'str',
            name:       'str',
            status:     'str',
            exit_status: 'num',
            submit:     'date',
            start:      'date',
            complete:   'date',
            wall_time:  'time',
            run_time:   'time',
            state:      'str',
            '%cpu':     'num',
            '%mem':     'num',
            rss:        'mem',
            vmem:       'mem',
            'max_%cpu': 'num',
            'max_%mem': 'num',
            max_rss:    'mem',
            max_vmem:   'mem',
            rchar:      'mem',
            wchar:      'mem',
            syscr:      'mem',
            syscw:      'mem',
            read_bytes: 'mem',
            write_bytes: 'mem',
            vmpeak:     'mem',
            vmsize:     'mem',
            vmlck:      'mem',
            vmpin:      'mem',
            vmhwm:      'mem',
            vmrss:      'mem',
            vmdata:     'mem',
            vmstk:      'mem',
            vmexe:      'mem',
            vmlib:      'mem',
            vmpte:      'mem',
            vmswap:     'mem',
            threads:    'num'
    ]

    static Map<String,Closure> FORMATTER = [
            str: this.&fmtToString,
            num: this.&fmtToNumber,
            date: this.&fmtToDate,
            time: this.&fmtToTime,
            mem: this.&fmtToMem
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

    @PackageScope
    static def fmtToString(def value, String fmt) {
        value ? value.toString() : NA
    }

    @PackageScope
    static def fmtToNumber(def value, String fmt) {
        value == Integer.MAX_VALUE ? NA :value.toString()
    }

    @PackageScope
    static def fmtToDate(def value, String fmt) {
        if( !value ) return NA
        getDateFormat(fmt).format(new Date(value as long))
    }

    @PackageScope
    static def fmtToTime(def value, String fmt) {
        if( value == null ) return NA
        new Duration(value as long).toString()
    }

    @PackageScope
    static def fmtToMem( def value, String fmt) {
        if( value == null ) return NA

        String str = value.toString()
        if( str.isLong() ) {
            str = new MemoryUnit(str.toLong()).toString()
            str = str.replaceAll(/,/,'')
        }

        return str
    }

//    @PackageScope
//    static def fmtToPerc( def value, String fmt ) {
//        try {
//            def num = value instanceof Number ? value.toFloat() : value.toString().toFloat()
//            return String.format('%.2f', num)
//        }
//        catch( Exception e ) {
//            log.debug "Not a valid percentual value: '$value'"
//            return '-'
//        }
//    }


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
        store.put(name, value)

        // when then 'complete' timestamp is set, calculates
        // the ''wall_time' and 'run_time' fields
        if( name == 'complete' && value instanceof Long ) {
            if( store.submit )
                store.wall_time = (value as long) - (store.submit as long)

            if( store.start )
                store.run_time = (value as long) - (store.start as long)
        }

        // vmpeak: Peak virtual memory size
        // this is a synonym of 'max_vmem' field
        else if( name == 'vmpeak' ) {
            store.put('max_vmem', value)
        }

        // Peak resident set size ("high water mark")
        // This is a synonym of 'max_rss' field
        else if( name == 'vmhwm' ) {
            store.put('max_rss', value)
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
    def get( String name, String converter ) {
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
            throw  new IllegalArgumentException("Not a valid trace formatter for field: '$name' with type: '$type'")

        formatter.call(val,sFormat)
    }

    def getTaskId() { get('task_id') }

    /**
     * Render the specified list of fields to a single string value
     *
     * @param fields The list of fields to be rendered, each entry can optionally specify a
     *      a format string separating it from the name by a colon character e.g. name:format
     * @param delim A delimiter that separates fields entries in the final string
     * @return The final string containing all field values
     */
    String render( List<String> fields, String delim ) {
        def result = new ArrayList(fields.size())
        for( def item : fields ) {
            String name
            String format
            int p = item.indexOf(':')
            if( p == -1 ) {
                name = item
                format = null
            }
            else {
                name = item.substring(0,p)
                format = item.substring(p+1)
            }

            result << ( get(name, format) ?: NA )
        }

        return result.join(delim)
    }


    String toString() {
        "${this.class.simpleName}${store}"
    }

    boolean equals(o) {
        if (this.is(o)) return true
        if (getClass() != o.class) return false

        TraceRecord that = (TraceRecord) o

        if (store != that.store) return false

        return true
    }

    int hashCode() {
        return store.hashCode()
    }
}