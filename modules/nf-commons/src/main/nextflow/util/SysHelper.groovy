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

package nextflow.util
import java.lang.management.ManagementFactory
import java.text.SimpleDateFormat

import com.sun.management.OperatingSystemMXBean
import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.file.FileHelper
/**
 * System helper methods
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class SysHelper {

    private static String DATE_FORMAT = 'dd-MMM-yyyy HH:mm'

    /**
     * Given a timestamp as epoch time convert to a string representation
     * using the {@link #DATE_FORMAT}
     *
     * @param dateInMillis
     *      The date as number of milliseconds from 1 Jan 1970
     * @param tz
     *      The {@link TimeZone} to use to format the date, or {@code null} to use the default time-zone
     * @return
     *      The formatted date string
     */
    static String fmtDate(long dateInMillis, TimeZone tz=null) {
        fmtDate(new Date(dateInMillis), tz)
    }

    /**
     * Given a {@link Date} object convert to a string representation
     * according the {@link #DATE_FORMAT}
     *
     * @param date
     *      The date to render as a string
     * @param tz
     *      The {@link TimeZone} to use to format the date, or {@code null} to use the default time-zone
     * @return
     *      The formatted date string
     */
    static String fmtDate(Date date, TimeZone tz=null) {
        def formatter=new SimpleDateFormat(DATE_FORMAT)
        if(tz) formatter.setTimeZone(tz)
        formatter.format(date)
    }

    /**
     * Read the system uptime as returned by the {@code uptime} Linux tool
     *
     * @return The {@code uptime} stdout text
     * @throws IllegalStateException If the {@code uptime} command return a non-zero exit status
     */
    static String getUptimeText() throws IllegalStateException {

        def proc = new ProcessBuilder('uptime').start()
        def status = proc.waitFor()
        if( status == 0 ) {
            return proc.text.trim()
        }

        throw new IllegalStateException("Unable to run system `uptime` command: exit=$status")
    }

    /**
     * The system uptime
     *
     * @return A {@link Duration} object representing the uptime of the system
     */
    static Duration getUptimeDuration() {
        def text = getUptimeText()
        def result = parseUptimeText(text)
        log.trace "Uptime $result -- parsed text: $text"
        return result
    }

    static long getBootTimeMillis() {
        System.currentTimeMillis() - getUptimeDuration().toMillis()
    }

    static String getBootTimeString() {
        new SimpleDateFormat(DATE_FORMAT).format(new Date(getBootTimeMillis()))
    }

    /**
     * Parse `uptime` command line tool output
     *
     * @param str The text stdout produced by the {@code uptime} system tool
     * @return
     * @throws IllegalArgumentException
     */
    @PackageScope
    static Duration parseUptimeText( String str ) throws IllegalArgumentException {
        try {
            def p = str.indexOf(' up ')
            def items = str.substring(p+4).tokenize(',')
            def uptime = items[0].trim().replace('mins','min').replace('hrs','hour')
            if( uptime.contains(':') ) {
                uptime = reformatHourAndMin(uptime)
            }
            else if( items.size()>1 && items[1].contains(":")) {
                uptime += " ${reformatHourAndMin(items[1])}"
            }
            return Duration.of(uptime)

        }
        catch( Exception e ) {
            throw new IllegalArgumentException("Not a valid uptime text: `$str`", e)
        }
    }

    private static String reformatHourAndMin(String text) {
        def items = text.tokenize(':')
        "${items[0]}h ${items[1]}m"
    }

    /**
     * @return The actual free space in the node local storage
     */
    static MemoryUnit getAvailDisk() {
        final free = FileHelper.getLocalTempPath().toFile().getFreeSpace()
        new MemoryUnit(free)
    }

    @Memoized
    static private OperatingSystemMXBean getSystemMXBean() {
        (OperatingSystemMXBean) ManagementFactory.getOperatingSystemMXBean()
    }

    /**
     * @return The number of CPUs available
     */
    static int getAvailCpus() {
        final int result = getSystemMXBean().getAvailableProcessors()
        return result
    }

    /**
     * @return The total system memory available
     */
    static MemoryUnit getAvailMemory() {
        new MemoryUnit(getSystemMXBean().getTotalPhysicalMemorySize())
    }

    /**
     * @return The system host name
     */
    static String getHostName() {
        System.getenv('HOSTNAME') ?: InetAddress.getLocalHost().getHostName()
    }


}
