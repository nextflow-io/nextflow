package nextflow.util

import groovy.transform.CompileDynamic
import groovy.transform.CompileStatic

import org.apache.commons.lang3.StringUtils

@CompileStatic
class BasicBolts {
    /**
     * Find all the best matches for the given example string in a list of values
     *
     * @param sample The example string -- cannot be empty
     * @param options A list of string
     * @return The list of options that best matches to the specified example -- return an empty list if none match
     */
    @CompileDynamic
    static List<String> closest(Collection<String> options, String sample ) {
        assert sample

        if( !options )
            return Collections.emptyList()

        // Otherwise look for the most similar
        def diffs = [:]
        options.each {
            diffs[it] = StringUtils.getLevenshteinDistance(sample, it)
        }

        // sort the Levenshtein Distance and get the fist entry
        def sorted = diffs.sort { it.value }
        def nearest = sorted.find()
        def min = nearest.value
        def len = sample.length()

        def threshold = len<=3 ? 1 : ( len > 10 ? 5 : Math.floor(len/2))

        def result
        if( min <= threshold ) {
            result = sorted.findAll { it.value==min } .collect { it.key }
        }
        else {
            result = []
        }

        return result
    }
}
