package nextflow.util

import groovy.transform.EqualsAndHashCode

/**
 * DnaNexus file handler
 *
 */

@EqualsAndHashCode(includes = ['id','name'])
class DxFile {

    String id

    String name


    String toString() {
        name
    }



}
