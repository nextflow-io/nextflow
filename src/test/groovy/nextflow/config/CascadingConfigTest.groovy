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
