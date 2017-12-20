package nextflow.mail

import java.nio.file.Path

import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString

/**
 * Model a mail attachment
 */
@ToString(includeNames = true)
@EqualsAndHashCode
class Attachment {
    /**
     * The attachment file
     */
    private File file

    /**
     * The attachment as classpath resource
     */
    private String resource

    /**
     * Attachment content parameters
     */
    private Map params

    static Attachment resource(Map params, String path) {
        def result = new Attachment()
        result.resource = path
        result.params = params != null ? params : [:]
        return result
    }

    protected Attachment() {}

    Attachment( attach ) {
        this([:], attach)
    }

    Attachment( Map params, attach ) {
        if( attach instanceof File ) {
            this.file = attach
        }
        else if( attach instanceof Path ) {
            this.file = attach.toFile()
        }
        else if( attach instanceof String || attach instanceof GString ) {
            this.file = new File(attach.toString())
        }
        else if( attach != null )
            throw new IllegalArgumentException("Invalid attachment argument: $attach [${attach.getClass()}]")

        this.params = params != null ? params : [:]
    }

    File getFile() { file }

    String getResource() { resource }

    String getContentId() { params.contentId }

    String getDescription() { params.description }

    String getDisposition() { params.disposition }

    String getFileName() {
        if( params.fileName )
            return params.fileName

        if( file )
            return file.getName()

        if( resource ) {
            def str = resource
            def p = str.indexOf(':')
            if( p!=-1 )
                str = resource.substring(p+1)
            p = str.lastIndexOf('/')
            if( p!=-1 )
                str = str.substring(p+1)
            return str
        }

        return null
    }
}