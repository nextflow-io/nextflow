/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
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

package nextflow.k8s

import nextflow.cli.CliOptions
import nextflow.cli.CmdRun
import nextflow.cli.Launcher
import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class K8sDriverLauncherTest extends Specification {

    @Unroll
    def 'should get cmd cli' () {

        given:
        def l = new K8sDriverLauncher(cmd: cmd, pipelineName: 'foo')

        when:
        cmd.launcher = new Launcher(options: new CliOptions())
        then:
        l.getLaunchCli() == expected

        where:
        cmd                                         | expected
        new CmdRun()                                | 'nextflow run foo'
        new CmdRun(cacheable: false)                | 'nextflow run foo -cache false'
        new CmdRun(resume: true)                    | 'nextflow run foo -resume true'
        new CmdRun(poolSize: 10)                    | 'nextflow run foo -ps 10'
        new CmdRun(pollInterval: 5)                 | 'nextflow run foo -pi 5'
        new CmdRun(queueSize: 9)                    | 'nextflow run foo -qs 9'
        new CmdRun(revision: 'xyz')                 | 'nextflow run foo -r xyz'
        new CmdRun(latest: true)                    | 'nextflow run foo -latest true'
        new CmdRun(withTrace: true)                 | 'nextflow run foo -with-trace true'
        new CmdRun(withTimeline: true)              | 'nextflow run foo -with-timeline true'
        new CmdRun(withDag: true)                   | 'nextflow run foo -with-dag true'
        new CmdRun(profile: 'ciao')                 | 'nextflow run foo -profile ciao'
        new CmdRun(dumpHashes: true)                | 'nextflow run foo -dump-hashes true'
        new CmdRun(dumpChannels: 'lala')            | 'nextflow run foo -dump-channels lala'
        new CmdRun(env: [XX:'hello', YY: 'world'])  | 'nextflow run foo -e.XX hello -e.YY world'
        new CmdRun(process: [mem: '100',cpus:'2'])  | 'nextflow run foo -process.mem 100 -process.cpus 2'
        new CmdRun(params: [alpha:'x', beta:'y'])   | 'nextflow run foo --alpha x --beta y'
    }

    @Unroll
    def 'should get pod name' () {

        given:
        def l = new K8sDriverLauncher(runName: name)

        expect:
        l.getPodName() == expect
        where:
        name        | expect
        'foo'       | 'nf-run-foo'
        'foo_bar'   | 'nf-run-foo-bar'
    }
}
