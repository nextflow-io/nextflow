/*
 * Copyright (c) 2020-2021. Seqera Labs, S.L.
 *
 * All Rights reserved
 *
 */

package nextflow.ui.console

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ConsoleRunnerTest extends Specification{

    def 'should load nextflow icon' () {
        expect:
        ConsoleRunner.getResource("/nextflow-icon.png")
    }

}
