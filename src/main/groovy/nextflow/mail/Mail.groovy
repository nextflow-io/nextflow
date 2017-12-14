package nextflow.mail

import java.nio.file.Path

import groovy.transform.EqualsAndHashCode

/**
 * Helper class modeling mail parameters
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@EqualsAndHashCode
class Mail {

    String from

    String to

    String cc

    String bcc

    String subject

    String charset

    String body

    String text

    List<File> attachments

    String type

    /**
     * Creates a {@link Mail} object given a {@link Mail} object
     *
     * @param params
     * @return A mail object representing the message to send
     */
    static Mail of(Map params) {
        def result = new Mail()

        if( params.from )
            result.from(params.from.toString())

        if( params.to )
            result.to(params.to.toString())

        if( params.cc )
            result.cc(params.cc.toString())

        if( params.bcc )
            result.bcc(params.bcc.toString())

        if( params.subject )
            result.subject(params.subject.toString())

        if( params.charset )
            result.charset(params.charset.toString())

        if( params.type )
            result.type(params.type.toString())

        if( params.body )
            result.body(params.body)

        if( params.text )
            result.text(params.text)

        if( params.attach )
            result.attach(params.attach)

        return result
    }

    /**
     * Mail sender (addresses must follow RFC822 syntax)
     * Multiple addresses can be separated by a comma
     */
    void from( String address ) {
        this.from = address
    }

    /**
     * Mail TO recipients (addresses must follow RFC822 syntax).
     * Multiple addresses can be separated by a comma.
     */
    void to( String address ) {
        this.to = address
    }

    /**
     * Mail CC recipients (addresses must follow RFC822 syntax).
     * Multiple addresses can be separated by a comma.
     */
    void cc( String address ) {
        this.cc = address
    }

    /**
     * Mail BCC recipients (addresses must follow RFC822 syntax).
     * Multiple addresses can be separated by a comma.
     */
    void bcc( String address ) {
        this.bcc = address
    }

    /**
     * @param subject The email subject
     */
    void subject( String subject ) {
        this.subject = subject
    }

    private String stringify( value ) {
        if( value instanceof File )
            return value.text
        if( value instanceof Path )
            return value.text
        if( value instanceof CharSequence )
            return value.toString()
        if( value != null )
            throw new IllegalArgumentException("Not a valid mail body argument [${value.getClass().getName()}]: $value")
        return null
    }

    /**
     * @param str The email content
     */
    void body( value ) {
        this.body = stringify(value)
    }

    /**
     * Plain text mail content
     * @param text The mail text content
     */
    void text( value ) {
        this.text = stringify(value)
    }

    /**
     * Mail content mime-type
     */
    void type( String mime )  {
        this.type = mime
    }

    /**
     * @param charset The mail content charset
     */
    void charset( String charset ) {
        this.charset = charset
    }

    /**
     * Add an email attachment
     *
     * @param item A attachment file path either as {@link File}, {@code Path} or {@link String}
     */
    void attach( item ) {
        if( this.attachments == null )
            this.attachments = []

        if( item instanceof Collection ) {
            item.each { attach0(it) }
        }
        else if( item instanceof Object[] ) {
            item.each { attach0(it) }
        }
        else if( item ) {
            attach0(item)
        }
    }

    /**
     * Add a single item to the list of file attachments
     * @param attach either a {@link File}, {@link java.nio.file.Path} or a string file path
     */
    private void attach0( attach ) {
        if( attach instanceof File ) {
            this.attachments.add(attach)
        }
        else if( attach instanceof Path ) {
            this.attachments.add(attach.toFile())
        }
        else if( attach instanceof String || attach instanceof GString ) {
            this.attachments.add(new File(attach.toString()))
        }
        else if( attach != null )
            throw new IllegalArgumentException("Invalid attachment argument: $attach [${attach.getClass()}]")
    }

}
