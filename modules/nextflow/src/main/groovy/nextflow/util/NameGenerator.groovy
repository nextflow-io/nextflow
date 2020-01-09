/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

package nextflow.util

import groovy.transform.CompileStatic
import groovy.transform.PackageScope

/**
 * Generate a random mnemonic name
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class NameGenerator {

    @PackageScope
    static List<String> ADJECTIVES = [

            "admiring",
            "adoring",
            "agitated",
            "amazing",
            "angry",
            "astonishing",
            "awesome",
            "backstabbing",
            "berserk",
            "big",
            "boring",
            "chaotic",
            "cheeky",
            "cheesy",
            "clever",
            "compassionate",
            "condescending",
            "confident",
            "cranky",
            "crazy",
            "curious",
            "deadly",
            "desperate",
            "determined",
            "distracted",
            "distraught",
            "disturbed",
            "dreamy",
            "drunk",
            "ecstatic",
            "elated",
            "elegant",
            "evil",
            "exotic",
            "extravagant",
            "fabulous",
            "fervent",
            "festering",
            "focused",
            "friendly",
            "furious",
            "gigantic",
            "gloomy",
            "golden",
            "goofy",
            "grave",
            "happy",
            "high",
            "hopeful",
            "hungry",
            "infallible",
            "insane",
            "intergalactic",
            "irreverent",
            "jolly",
            "jovial",
            "kickass",
            "lethal",
            "lonely",
            "loving",
            "mad",
            "magical",
            "maniac",
            "marvelous",
            "mighty",
            "modest",
            "nasty",
            "naughty",
            "nauseous",
            "nice",
            "nostalgic",
            "peaceful",
            "pedantic",
            "pensive",
            "prickly",
            "reverent",
            "ridiculous",
            "romantic",
            "sad",
            "scruffy",
            "serene",
            "sharp",
            "shrivelled",
            "sick",
            "silly",
            "sleepy",
            "small",
            "soggy",
            "special",
            "spontaneous",
            "stoic",
            "stupefied",
            "suspicious",
            "tender",
            "thirsty",
            "tiny",
            "trusting",
            "voluminous",
            "wise",
            "zen"

    ]

    @PackageScope
    static final List<String> NAMES = [
//  All TaylorLab Members
            "alex",
            "allison",
            "evan",
            "mark.d",
            "nicole",
            "noah",
            "pamela",
            "sam",
            "shweta",
            "tripti",
            "vasilisa",
            "yixiao",
            "barry",
            "chai",
            "craig",
            "ezra",
            "helen",
            "mark.z",
            "phillip",
            "tambu

    ]

    private static Random RND = new Random()

    /**
    * @return A random generated name string
    */
    static String next() {
        def a = ADJECTIVES[ RND.nextInt(ADJECTIVES.size()) ]
        def b = NAMES[ RND.nextInt(NAMES.size()) ]
        return "${a}_${b}"
    }

    static String next(String... skip) {
        next(skip as List<String>)
    }

    static String next(Collection<String> skip) {
        while( true ) {
            final result = next()
            if( !skip.contains(result) )
                return result
        }
    }

}
