/*
 * Copyright 2013-2024, Seqera Labs
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
import java.nio.charset.Charset

import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CharsetHelperTest extends Specification {

    def testGetCharSet() {

        expect:
        CharsetHelper.getCharset('x') == Charset.defaultCharset()
        CharsetHelper.getCharset([x:1]) == Charset.defaultCharset()
        CharsetHelper.getCharset('iso-8859-1') == Charset.forName('iso-8859-1')
        CharsetHelper.getCharset(charset:'iso-8859-1') == Charset.forName('iso-8859-1')
        CharsetHelper.getCharset(xx:'iso-8859-1') == Charset.defaultCharset()

    }

}
