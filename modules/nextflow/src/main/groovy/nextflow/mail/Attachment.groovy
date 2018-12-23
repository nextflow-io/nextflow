/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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