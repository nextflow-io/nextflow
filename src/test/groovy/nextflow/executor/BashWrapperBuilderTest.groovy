package nextflow.executor

import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class BashWrapperBuilderTest extends Specification {

    def 'test changeToScratchDir' () {

        setup:
        def builder = [:] as BashWrapperBuilder

        expect:
        builder.changeToScratchDirectory() == null

        when:
        builder.scratch = true
        then:
        builder.changeToScratchDirectory() == 'NF_SCRATCH=${TMPDIR:-`mktemp -d`} && cd $NF_SCRATCH'

        when:
        builder.scratch = '$SOME_DIR'
        then:
        builder.changeToScratchDirectory() == 'NF_SCRATCH=${SOME_DIR:-`mktemp -d`} && cd $NF_SCRATCH'

        when:
        builder.scratch = '/my/temp'
        then:
        builder.changeToScratchDirectory() == 'NF_SCRATCH=$(mktemp -d -p /my/temp) && cd $NF_SCRATCH'

        when:
        builder.scratch = '/my/temp'
        then:
        builder.changeToScratchDirectory() == 'NF_SCRATCH=$(mktemp -d -p /my/temp) && cd $NF_SCRATCH'

    }



}
