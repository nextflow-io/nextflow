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

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class IdentityTest extends Specification{

    def 'should validate hash code' () {

        expect:
        new HexIdentity(1).hashCode() == new HexIdentity(1).hashCode()
        new HexIdentity(100).hashCode() == new HexIdentity(100).hashCode()
        new HexIdentity(100).hashCode() != new HexIdentity(101).hashCode()
    }


    def "should validate equals methods" () {

        expect:
        new HexIdentity(1) == new HexIdentity(1)
        new HexIdentity(100) == new HexIdentity(100)
        new HexIdentity(100) != new HexIdentity(101)

    }

    def "should validate toString method" () {

        expect:
        new HexIdentity(1).toString() == "0x1"
        new HexIdentity(10).toString() == "0xa"
        new HexIdentity(100).toString() == "0x64"

    }

}
