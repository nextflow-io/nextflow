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

package nextflow.config

import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CascadingConfigTest extends Specification {

    static class FakeConfig extends CascadingConfig<String,Object> {

        FakeConfig(Map<String,Integer> cfg, CascadingConfig<String,Object> fallback)  {
            super(cfg,fallback)
        }

        @ConfigField
        def getFoo() {
            getAttribute('foo')
        }

        @ConfigField
        def getBar() {
            getAttribute('bar')
        }


        @ConfigField('user')
        def getUserName() {
            getAttribute('user')
        }

    }


    def 'should retrieve values' () {

        given:
        def cfg = new FakeConfig([foo:1,bar:2], new FakeConfig([bar:3, user: 5], null))
        expect:
        cfg.getFoo() == 1
        cfg.getBar() == 2
        cfg.getUserName() == 5

        when:
        cfg.getAttribute('unknown')
        then:
        thrown(IllegalArgumentException)

        expect:
        cfg.validFields() == ['foo','bar','user'] as Set
    }

    def 'should return null attribute when missing' () {

        given:
        def fallback = new CascadingConfig() { }

        when:
        def cfg = new FakeConfig([foo:1,bar:2], fallback)
        then:
        cfg.getAttribute('foo') == 1
        cfg.getAttribute('bar') == 2
        cfg.getAttribute('user') == null


    }


}
