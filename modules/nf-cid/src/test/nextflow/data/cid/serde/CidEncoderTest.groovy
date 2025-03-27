package nextflow.data.cid.serde

import nextflow.data.cid.model.Checksum
import nextflow.data.cid.model.Output
import nextflow.data.cid.model.TaskOutput
import spock.lang.Specification

class CidEncoderTest extends Specification{

    def 'should encode and decode Outputs'(){
        given:
            def encoder = new CidEncoder()
        and:
            def output = new TaskOutput("/path/to/file", new Checksum("hash_value", "hash_algorithm", "standard"), "cid://source", 1234)

        when:
            def encoded = encoder.encode(output)
            def object = encoder.decode(encoded)

        then:
            object instanceof Output
            output.path == "/path/to/file"
            output.checksum instanceof Checksum
            output.checksum.value == "hash_value"
            output.checksum.algorithm == "hash_algorithm"
            output.checksum.mode == "standard"
            output.source == "cid://source"
            output.size == 1234

    }
}
