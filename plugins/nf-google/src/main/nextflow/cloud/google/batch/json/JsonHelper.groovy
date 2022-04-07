/*
 * Copyright (c) 2020-2021. Seqera Labs, S.L.
 *
 * All Rights reserved
 *
 */

package nextflow.cloud.google.batch.json

import groovy.json.JsonGenerator
import groovy.json.JsonOutput
import groovy.transform.CompileStatic

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class JsonHelper {

    private static final JsonGenerator generator

    static {
        generator = new JsonGenerator.Options()
                .excludeNulls()
                .build()
    }

    static String toJson(Object o, boolean pretty=false) {
        final json = generator.toJson(o)
        return pretty ? JsonOutput.prettyPrint(json) : json
    }

}
