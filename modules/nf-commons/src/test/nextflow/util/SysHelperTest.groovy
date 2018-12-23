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

package nextflow.util

import java.lang.management.ManagementFactory

import com.sun.management.OperatingSystemMXBean
import nextflow.file.FileHelper
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class SysHelperTest extends Specification {


    def 'should parse uptime to duration' () {
        expect:
        SysHelper.parseUptimeText(str) == duration
        where:
        str | duration
        '20:29  up 2 days,  8:33, 4 users, load averages: 4.41 3.05 2.65'           | Duration.of('2d 8h 33min')
        ' 14:39:19 up 430 days,  4:16, 28 users,  load average: 0.42, 0.19, 0.05'   | Duration.of('430d 4h 16min')
        ' 20:28:46 up 129 days,  3:16,  2 users,  load average: 0.00, 0.01, 0.05'   | Duration.of('129d 3h 16min')
        ' 19:22:07 up 4 min,  1 user,  load average: 0.01, 0.08, 0.05'              | Duration.of('4min')
        ' 9:21  up 35 mins, 5 users, load averages: 3.57 2.48 2.00'                 | Duration.of('35min')
        '10:04  up  1:18, 5 users, load averages: 3.94 3.64 3.68'                   | Duration.of('1h 18m')
        ' 08:05:44 up 0 min,  1 user,  load average: 1.10, 0.27, 0.09'              | Duration.of('0m')
        '13:24  up 17 hrs, 3 users, load averages: 5.00 4.33 3.67'                  | Duration.of('17 hour')
        '13:27  up 17:03, 3 users, load averages: 4.27 4.11 3.68'                   | Duration.of('17 h 3 min')
        '14:31  up 1 day, 18:08, 7 users, load averages: 2.47 3.10 3.57'            | Duration.of('1d 18h 8m')

    }

    def 'should read uptime' () {

        expect:
        SysHelper.getUptimeDuration().toMillis() > 0

    }

    def 'should get boot time' () {
        given:
        def bootTime = SysHelper.getBootTimeMillis()
        println "boot: ${SysHelper.getBootTimeString()} [${SysHelper.getUptimeDuration()}]"
        expect:
        bootTime > 0
        bootTime < System.currentTimeMillis()
    }

    def 'should validate sys metrics' () {
        given:
        def bean = (OperatingSystemMXBean) ManagementFactory.getOperatingSystemMXBean()
        expect:
        SysHelper.availCpus == Runtime.runtime.availableProcessors()
        SysHelper.availMemory == new MemoryUnit(bean.getTotalPhysicalMemorySize())
        SysHelper.availDisk.giga == new MemoryUnit(FileHelper.getLocalTempPath().toFile().getFreeSpace()).giga
        SysHelper.hostName == (System.getenv('HOSTNAME') ?: InetAddress.getLocalHost().getHostName())
    }

    def 'should format date string' () {
        expect:
        SysHelper.fmtDate(1470901220000, TimeZone.getTimeZone('GMT+2')) == '11-Aug-2016 09:40'
    }

}
