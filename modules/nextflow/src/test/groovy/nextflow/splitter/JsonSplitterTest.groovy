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
import com.google.gson.stream.JsonReader;

/**
 * @author Pierre Lindenbaum PhD Institut-du-Thorax Nantes France.
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
                "inserted_sequence": ["G",["C","A"]]
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

    def testSplitArrayLimit() {
        when:
        def items = new JsonSplitter(limit:2).target(json_array).list()
        
        then:
        items.size() == 2
        items[0].get("beaconId").equals("se.nbis.swefreq");
    }


    def testSplitObject() {
        when:
        def items = new JsonSplitter().target(json_object).list()

        then:
        items.size() == 4
        items[0] instanceof Map
        items[0].get(JsonSplitter.OBJECT_KEY).equals("seq_id");
        items[0].get(JsonSplitter.OBJECT_VALUE).equals("NC_000008.9");
        items[1].get(JsonSplitter.OBJECT_KEY).equals("position");
        items[1].get(JsonSplitter.OBJECT_VALUE) == 19864003
    }

    def testSplitObjectLimit() {
        when:
        def items = new JsonSplitter(limit:1).target(json_object).list()

        then:
        items.size() == 1
        items[0].get(JsonSplitter.OBJECT_KEY).equals("seq_id");
        items[0].get(JsonSplitter.OBJECT_VALUE).equals("NC_000008.9");
    }


    def testJsonPathArray1() {
        when:
        def items = new JsonSplitter("jsonPath":"[0].alleleRequest").
        target(json_array).list()

        then:
        items.size() == 6
        items[0].get(JsonSplitter.OBJECT_KEY).equals("referenceName");
        items[0].get(JsonSplitter.OBJECT_VALUE).equals("MT");
    }


    def testJsonPathObject1() {
        when:
        def items = new JsonSplitter("jsonPath":"inserted_sequence[1]").
        target(json_object).list()

        then:
        items.size() == 2
        items[0].equals("C")
        items[1].equals("A")
    }


    def testFromJson() {
        when:
        StringReader r = new StringReader(
                "[true,false,null,1,10000000000,1E1000,1.0,[],{}]"
                );
        JsonReader jr = new JsonReader(r);
        def items = JsonSplitter.fromJson(jr);
        jr.close();
        r.close();
        
        then:
        items instanceof List
        items[0] == true
        items[1] == false
        items[2] == null
        items[3] instanceof Integer
        items[4] instanceof Long
        items[5] instanceof String
        items[6] instanceof Double
        items[7] instanceof List
        items[8] instanceof Map
    }
}
