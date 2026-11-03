/*
 * Copyright 2013-2026, Seqera Labs
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

package nextflow.cli

import nextflow.exception.AbortOperationException
import spock.lang.Specification

class CmdRemovedTest extends Specification {

    def 'should report that the command was removed'() {
        when:
        cmd.run()
        then:
        def e = thrown(AbortOperationException)
        e.message.startsWith(MESSAGE)

        when:
        cmd.usage()
        then:
        thrown(AbortOperationException)

        when:
        cmd.usage([])
        then:
        thrown(AbortOperationException)

        where:
        cmd                 | MESSAGE
        new CmdConsole()    | 'The `console` command has been removed'
        new CmdKubeRun()    | 'The `kuberun` command has been removed -- use `nextflow plugin nf-kuberun:run`'
    }

    def 'should parse leftover arguments of removed commands'() {
        when:
        def launcher = new Launcher().parseMainArgs(args as String[])
        then:
        launcher.command.class == CLASS

        where:
        args                                                | CLASS
        ['console']                                         | CmdConsole
        ['console', 'foo.nf']                               | CmdConsole
        ['console', '-foo']                                 | CmdConsole
        ['console', '-foo', 'bar.nf']                       | CmdConsole
        ['kuberun']                                         | CmdKubeRun
        ['kuberun', 'nextflow-io/hello']                    | CmdKubeRun
        ['kuberun', '-n', 'ns', 'nextflow-io/hello']        | CmdKubeRun
        ['kuberun', 'nextflow-io/hello', '-v', 'pvc:/mnt']  | CmdKubeRun
    }
}
