/*
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

package nextflow.splitter

import spock.lang.Specification
import spock.lang.Unroll

import java.nio.charset.Charset

/**
 *
 * @author Pierre Lindenbaum PhD
 */
class JsonSplitterTest extends Specification {

    def json_array = '''
        	[
		    {
			"beaconId": "se.nbis.swefreq",
			"apiVersion": 1.0,
			"exists": false,
			"alleleRequest": {
			    "referenceName": "MT",
			    "referenceBases": "T",
			    "assemblyId": "GRCh38",
			    "includeDatasetResponses": "NONE",
			    "alternateBases": "C",
			    "start": 9
			},
			"datasetAlleleResponses": [],
			"beaconHandover": []
		    },
		    {
			"beaconId": "fi.rahtiapp.beaconpy-elixirbeacon",
			"apiVersion": 1,
			"exists": true,
			"alleleRequest": {
			    "referenceName": "MT",
			    "referenceBases": "T",
			    "assemblyId": "GRCh38",
			    "includeDatasetResponses": "NONE",
			    "alternateBases": "C",
			    "start": 9
			},
			"datasetAlleleResponses": []
		    },
		    {
			"beaconId": "fi.rahtiapp.staging-elixirbeacon",
			"apiVersion": "1.1.0",
			"exists": true,
			"alleleRequest": {
			    "referenceName": "MT",
			    "referenceBases": "T",
			    "assemblyId": "GRCh38",
			    "includeDatasetResponses": "NONE",
			    "alternateBases": "C",
			    "start": 9
			},
			"datasetAlleleResponses": []
		    }
		]
        '''

    def json_object = '''
    	    {
                "seq_id": "NC_000008.9",
                "position": 19864003,
                "deleted_sequence": "C",
                "inserted_sequence": "G"
            }
    	'''


    def testSplitArray() {

        when:
		def items = new JsonSplitter().target(json_array).list()
        then:
		items.size() == 3
		items[0] instanceof Map
		items[0].get("beaconId").equals("se.nbis.swefreq");
		items[0].get("exists") instanceof Boolean
		items[0].get("apiVersion") instanceof Double
		items[1].get("apiVersion") instanceof Integer

    }
    
    def testSplitObject() {

        when:
        	def items = new JsonSplitter().target(json_object).list()

        then:
		items.size() == 4
		items[0] instanceof Map
		items[0].get("seq_id").equals("NC_000008.9");
		items[1].get("position") == 19864003
    }
    
    def testNotAJsonFile() {

        when:
        	def items = new JsonSplitter().target("[not a json string").list()
        then:
      		 thrown(com.google.gson.stream.MalformedJsonException)
    }
    
}
