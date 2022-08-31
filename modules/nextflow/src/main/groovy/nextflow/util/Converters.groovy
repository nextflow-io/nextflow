package nextflow.util

import groovy.json.JsonBuilder
import groovy.json.JsonOutput
import org.codehaus.groovy.runtime.InvokerHelper

/**
 * @author : jorge <jorge.aguilera@seqera.io>
 *
 * Some generic converters methods
 */
class Converters {

    static def deepReplaceToString(root, replaceNullWith = "") {
        if (root instanceof List) {
            root.collect {
                if (it instanceof Map) {
                    deepReplaceToString(it, replaceNullWith)
                } else if (it instanceof List) {
                    deepReplaceToString(it, replaceNullWith)
                } else if (it == null) {
                    replaceNullWith
                } else {
                    it.toString()
                }
            }
        } else if (root instanceof Map) {
            root.each {
                if (it.value instanceof Map) {
                    deepReplaceToString(it.value, replaceNullWith)
                } else if (it.value instanceof List) {
                    it.value = deepReplaceToString(it.value, replaceNullWith)
                } else if (it.value == null) {
                    it.value = replaceNullWith
                } else {
                    it.value = it.value.toString()
                }
            }
        }
    }

    static String prettyPrint(input ){
        if (input instanceof Map) {
            def converted =  deepReplaceToString(input)
            return JsonOutput.prettyPrint(JsonOutput.toJson(converted))
        } else if (input instanceof List) {
            def converted = deepReplaceToString(input)
            return new JsonBuilder(converted).toPrettyString()
        } else {
            return InvokerHelper.inspect(input)
        }
    }

}
