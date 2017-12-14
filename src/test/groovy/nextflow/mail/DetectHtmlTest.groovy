package nextflow.mail

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DetectHtmlTest extends Specification {

    def 'should detect html content' () {

        expect:
        !DetectHtml.isHtml('Hello')
        !DetectHtml.isHtml('< hello')
        !DetectHtml.isHtml('< Hello>')
        DetectHtml.isHtml('<hello> xx </hello>')
        DetectHtml.isHtml('<br/>')

    }

}
