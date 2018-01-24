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