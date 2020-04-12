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

package nextflow.file

import groovy.json.JsonSlurper
import org.yaml.snakeyaml.Yaml
import spock.lang.Specification
import test.TestHelper

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class SlurperExTest extends Specification {

    def 'should parse json from path' () {
        given:
        def file = TestHelper.createInMemTempFile('foo.json')
        file.text = '''
        [
          {
            "patient_id": "ATX-TBL-001-GB-01-105",
            "region_id": "R1",
            "feature": "pass_vafqc_flag",
            "pass_flag": "TRUE"
          },
          {
            "patient_id": "ATX-TBL-001-GB-01-105",
            "region_id": "R1",
            "feature": "pass_stripy_flag",
            "pass_flag": "TRUE"
          } ]    
        '''
        
        when:
        def result = new JsonSlurper().parse(file)
        then:
        result.size() == 2
        and:
        result[0].patient_id == "ATX-TBL-001-GB-01-105"
        result[0].region_id == 'R1'
        result[0].feature == 'pass_vafqc_flag'
        and:
        result[1].patient_id == "ATX-TBL-001-GB-01-105"
        result[1].region_id == 'R1'
        result[1].feature == 'pass_stripy_flag'

    }

    def 'should parse yaml' () {
        given:
        def file = TestHelper.createInMemTempFile('foo.json')
        file.text = '''
-
  patient_id: ATX-TBL-001-GB-01-105
  region_id: R1
  feature: pass_vafqc_flag
  pass_flag: 'TRUE'
-
  patient_id: ATX-TBL-001-GB-01-105
  region_id: R1
  feature: pass_stripy_flag
  pass_flag: 'TRUE'
'''
        when:
        def result = new Yaml().load(file)
        
        then:
        result.size() == 2
        and:
        result[0].patient_id == "ATX-TBL-001-GB-01-105"
        result[0].region_id == 'R1'
        result[0].feature == 'pass_vafqc_flag'
        and:
        result[1].patient_id == "ATX-TBL-001-GB-01-105"
        result[1].region_id == 'R1'
        result[1].feature == 'pass_stripy_flag'
    }

}
