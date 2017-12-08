package nextflow.mail

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class MailParamsTest extends Specification {

    def 'should capture mail params' () {

        given:
        def closure = {
            from 'jim@dot.com'
            to 'paolo@dot.com'
            cc 'you@dot.com'
            bcc 'mrhide@dot.com'
            subject 'Hi there'
            content 'Hello world'
            attach 'foo.png'
            attach (['this.txt','that.txt'])

        }

        when:
        def params = new MailParams()
        params.with(closure)
        then:
        params.delegate.from == 'jim@dot.com'
        params.delegate.to == 'paolo@dot.com'
        params.delegate.cc == 'you@dot.com'
        params.delegate.bcc == 'mrhide@dot.com'
        params.delegate.subject ==  'Hi there'
        params.delegate.content == 'Hello world'
        params.delegate.attach == ['foo.png', 'this.txt','that.txt']


    }

}
