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

package nextflow.processor

import spock.lang.Specification
import test.TestHelper

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TaskHandlerTest extends Specification {

    def static final long KB = 1024
    
    def testParseTrace() {

        given:
        def traceText =  '''
        %cpu  %mem  rss  vmem  state
        2.0 1.0 968 2433364 R+
        4.0 3.0 1400 3000000
        rchar: 4786
        wchar: 12
        syscr: 11
        syscw: 1
        read_bytes: 0
        write_bytes: 0
        cancelled_write_bytes: 0
        Name:\tbash
        State:\tS (sleeping)
        Tgid:\t2515
        Pid:\t2515
        PPid:\t2513
        TracerPid:\t0
        Uid:\t1000\t1000\t1000\t1000
        Gid:\t1000\t1000\t1000\t1000
        FDSize:\t256
        Groups:\t1000
        VmPeak:\t   11084 kB
        VmSize:\t   11084 kB
        VmLck:\t       0 kB
        VmPin:\t       0 kB
        VmHWM:\t    1228 kB
        VmRSS:\t    1228 kB
        VmData:\t     128 kB
        VmStk:\t     136 kB
        VmExe:\t     900 kB
        VmLib:\t    2084 kB
        VmPTE:\t      40 kB
        VmSwap:\t       0 kB
        Threads:\t1
        SigQ:\t0/15951
        SigPnd:\t0000000000000000
        ShdPnd:\t0000000000000000
        SigBlk:\t0000000000010004
        SigIgn:\t0000000000000006
        SigCgt:\t0000000000010000
        CapInh:\t0000000000000000
        CapPrm:\t0000000000000000
        CapEff:\t0000000000000000
        CapBnd:\t0000001fffffffff
        Seccomp:\t0
        Cpus_allowed:\t1
        Cpus_allowed_list:\t0
        Mems_allowed:\t00000000,00000001
        Mems_allowed_list:\t0
        voluntary_ctxt_switches:\t32
        nonvoluntary_ctxt_switches:\t3
        '''
                .leftTrim()


        when:
        def handler = [:] as TaskHandler
        def trace = handler.parseTraceFile(traceText)
        then:
        trace.state == 'R+'
        trace.'%cpu' == '2.0'
        trace.'%mem' == '1.0'
        trace.rss == 968 * KB
        trace.vmem == 2433364 * KB
        trace.'max_%cpu' == '4.0'
        trace.'max_%mem' == '3.0'
        trace.max_rss == 1400 * KB
        trace.max_vmem == 3000000 * KB
        trace.vmpeak ==  11084 * KB
        trace.vmsize ==  11084 * KB
        trace.vmlck ==  0
        trace.vmpin ==  0
        trace.vmhwm ==  1228 * KB
        trace.vmrss ==  1228 * KB
        trace.vmdata == 128 * KB
        trace.vmstk ==  136 * KB
        trace.vmexe ==  900 * KB
        trace.vmlib ==  2084 * KB
        trace.vmpte ==  40 * KB
        trace.vmswap ==  0
        trace.threads ==  '1'
        trace.rchar == '4786'
        trace.wchar == '12'
        trace.syscr == '11'
        trace.syscw ==  '1'
        trace.read_bytes == '0'
        trace.write_bytes == '0'

    }


    def testGetTraceSmall() {

        given:
        def text = '''
        %cpu  %mem  rss  vmem  state
        1.0 2.0 968 2433364 S+
        3.0 4.0 2500 3503364
        '''.stripIndent().trim()

        def folder = TestHelper.createInMemTempDir()
        folder.resolve( TaskRun.CMD_TRACE ).text = text

        when:
        def handler = [:] as TaskHandler
        handler.task = new TaskRun(workDir: folder)
        def trace = handler.getTraceRecord()

        then:
        trace.get('state', null) == 'S+'
        trace.get('%cpu', null) == '1.0'
        trace.get('%mem', null) == '2.0'
        trace.get('rss', null ) == '968 KB'
        trace.get('vmem', null ) == '2.3 GB'
        trace.get('max_%cpu', null) == '3.0'
        trace.get('max_%mem', null) == '4.0'
        trace.get('max_rss', null ) == '2.4 MB'
        trace.get('max_vmem', null ) == '3.3 GB'
        trace.get('vmpeak', null ) == '-'
        trace.get('vmsize', null ) == '-'

    }

    def testGetTraceRecordFull() {

        given:
        def traceText =  '''
        %cpu  %mem  rss  vmem  state
        2.0 1.0 968 2433364 R+
        4.0 3.0 1400 3000000
        rchar: 4786
        wchar: 12
        syscr: 11
        syscw: 1
        read_bytes: 0
        write_bytes: 0
        cancelled_write_bytes: 0
        Name:\tbash
        State:\tS (sleeping)
        Tgid:\t2515
        Pid:\t2515
        PPid:\t2513
        TracerPid:\t0
        Uid:\t1000\t1000\t1000\t1000
        Gid:\t1000\t1000\t1000\t1000
        FDSize:\t256
        Groups:\t1000
        VmPeak:\t   11084 kB
        VmSize:\t   11084 kB
        VmLck:\t       0 kB
        VmPin:\t       0 kB
        VmHWM:\t    1228 kB
        VmRSS:\t    1228 kB
        VmData:\t     128 kB
        VmStk:\t     136 kB
        VmExe:\t     900 kB
        VmLib:\t    2084 kB
        VmPTE:\t      40 kB
        VmSwap:\t       0 kB
        Threads:\t1
        SigQ:\t0/15951
        SigPnd:\t0000000000000000
        ShdPnd:\t0000000000000000
        SigBlk:\t0000000000010004
        SigIgn:\t0000000000000006
        SigCgt:\t0000000000010000
        CapInh:\t0000000000000000
        CapPrm:\t0000000000000000
        CapEff:\t0000000000000000
        CapBnd:\t0000001fffffffff
        Seccomp:\t0
        Cpus_allowed:\t1
        Cpus_allowed_list:\t0
        Mems_allowed:\t00000000,00000001
        Mems_allowed_list:\t0
        voluntary_ctxt_switches:\t32
        nonvoluntary_ctxt_switches:\t3
        '''
                .leftTrim()

        def folder = TestHelper.createInMemTempDir()
        folder.resolve( TaskRun.CMD_TRACE ).text = traceText

        when:
        def handler = [:] as TaskHandler
        handler.task = new TaskRun(workDir: folder)
        def trace = handler.getTraceRecord()

        then:
        trace.state == 'R+'
        trace.'%cpu' == '2.0'
        trace.'%mem' == '1.0'
        trace.rss == 968 * KB
        trace.vmem == 2433364 * KB
        trace.'max_%cpu' == '4.0'
        trace.'max_%mem' == '3.0'
        trace.max_rss == 1228 * KB      // <-- note: this value is overridden by 'VmHWM'
        trace.max_vmem == 11084 * KB    // <-- note: this value is overridden by 'VmPeak'
        trace.vmpeak ==  11084 * KB
        trace.vmsize ==  11084 * KB
        trace.vmlck ==  0
        trace.vmpin ==  0
        trace.vmhwm ==  1228 * KB
        trace.vmrss ==  1228 * KB
        trace.vmdata == 128 * KB
        trace.vmstk ==  136 * KB
        trace.vmexe ==  900 * KB
        trace.vmlib ==  2084 * KB
        trace.vmpte ==  40 * KB
        trace.vmswap ==  0
        trace.threads ==  '1'
        trace.rchar == '4786'
        trace.wchar == '12'
        trace.syscr == '11'
        trace.syscw ==  '1'
        trace.read_bytes == '0'
        trace.write_bytes == '0'

        // check get method
        trace.get('state', null) == 'R+'
        trace.get('%cpu', null) == '2.0'
        trace.get('%mem', null) == '1.0'
        trace.get('rss', null ) == '968 KB'
        trace.get('vmem', null ) == '2.3 GB'

        trace.get('max_%cpu', null) == '4.0'
        trace.get('max_%mem', null) == '3.0'
        trace.get('max_rss', null ) == '1.2 MB'
        trace.get('max_vmem', null ) == '10.8 MB'

        trace.get('vmpeak', null ) == '10.8 MB' // (max_vmem)
        trace.get('vmhwm', null )  == '1.2 MB'  // (max_rss)
        trace.get('vmrss', null ) == '1.2 MB'   // (rss)
        trace.get('vmsize', null ) == '10.8 MB' // (vmem)

    }


}
