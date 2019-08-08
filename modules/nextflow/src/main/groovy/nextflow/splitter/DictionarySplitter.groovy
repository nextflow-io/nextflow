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

package nextflow.splitter

import java.nio.charset.Charset
import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.transform.InheritConstructors
import groovy.util.logging.Slf4j


import htsjdk.variant.utils.SAMSequenceDictionaryExtractor
import htsjdk.samtools.SAMSequenceDictionary
import htsjdk.samtools.SAMSequenceRecord

/**
 * Extracts SAM SequenceDictionary from SAM/BAM/CRAM/VCF/FASTA... files
 *
 * @author Pierre Lindenbaum PhD Institut du Thorax Nantes France.
 */
@Slf4j
@CompileStatic
@InheritConstructors
class DictionarySplitter extends AbstractTextSplitter {
   /* wich fields to fetch from a SAMSequenceRecord */
   private String ssrAttributes = SAMSequenceRecord.SEQUENCE_NAME_TAG+","+SAMSequenceRecord.SEQUENCE_LENGTH_TAG;
	
   @Override
   protected Reader newReader( Path path, Charset charset ) {
       List<String> lines = new ArrayList();
       final SAMSequenceDictionary dict = SAMSequenceDictionaryExtractor.extractDictionary(path);
       if(dict!=null) {
		def attributes = ssrAttributes.split("[,]");
		for(final SAMSequenceRecord ssr : dict.getSequences() ) {
			StringBuilder line = new StringBuilder();
			for(final String att: attributes) {
				String value;
				if(att.equals(SAMSequenceRecord.SEQUENCE_NAME_TAG)) value= ssr.getSequenceName();
				else if(att.equals(SAMSequenceRecord.SEQUENCE_LENGTH_TAG )) value= ssr.getSequenceLength();
				else if(att.equals(SAMSequenceRecord.ASSEMBLY_TAG  )) value= ssr.getAssembly();
				else value="";
				line.append(value).append("\t");
				}	
			line.append(path.toString()); 
			lines.add(line.toString());
			}
	}

        new BufferedReader(new StringReader(String.join("\n",lines)))
    }

    @Override
    protected Map<String,Object> validOptions() {
        def result = super.validOptions()
        result.attribute = [String]
        return result
    }

    @Override
    DictionarySplitter options(Map opts) {
        super.options(opts)
        if(opts.attribute) this.ssrAttributes = opts.attribute;
        return this
    }

    /**
     * A record is a text line
     *
     * @param reader The buffered reader
     * @return A line string or {@code null} when the end of the file is reached
     */
    @Override
    protected fetchRecord(BufferedReader reader) {
        def line = reader.readLine()
        return line==null? null : line.split("[\t]");
    }

}
