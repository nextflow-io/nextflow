package nextflow.mail

import groovy.transform.PackageScope

/**
 * Helper class modeling mail parameters
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class MailParams {

    private Map delegate = [:]

    @PackageScope Map getDelegate() { delegate }

    void from( address ) { delegate.from = address }

    void to( address ) { delegate.to = address }

    void cc( address ) { delegate.cc = address }

    void bcc( address ) { delegate.bcc = address }

    void subject( str ) { delegate.subject = str }

    void content( str ) { delegate.content = str }

    void attach( item ) {
        if( delegate.attach == null )
            delegate.attach = []

        if( item instanceof Collection ) {
            delegate.attach.addAll(item)
        }
        else if( item instanceof Object[] ) {
            delegate.attach.addAll(item)
        }
        else if( item ) {
            delegate.attach.add(item)
        }
    }
}
