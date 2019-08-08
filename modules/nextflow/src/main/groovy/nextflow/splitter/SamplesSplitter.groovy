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

import htsjdk.samtools.util.FileExtensions
import htsjdk.samtools.SamReaderFactory
import htsjdk.samtools.SAMFileHeader
import htsjdk.samtools.SAMReadGroupRecord
import htsjdk.samtools.ValidationStringency
import htsjdk.variant.vcf.VCFFileReader
import htsjdk.variant.vcf.VCFHeader


/**
 * Extracts Samples from SAM/BAM/CRAM/VCF files
 *
 * @author Pierre Lindenbaum PhD Institut du Thorax Nantes France.
 */
@Slf4j
@CompileStatic
@InheritConstructors
class SamplesSplitter extends AbstractTextSplitter {
   /* wich fields to fetch when using the BAM */
   private String rgAttributes = SAMReadGroupRecord.READ_GROUP_SAMPLE_TAG
	
   @Override
   protected Reader newReader( Path path, Charset charset ) {
        final String pathStr = path.toString();
        	List<String> lines = new ArrayList();
		if( pathStr.endsWith(FileExtensions.SAM) ||
			pathStr.endsWith(FileExtensions.BAM) 
			)
			{
			def attributes = this.rgAttributes.split("[,]");
			
			final Set<String> samples = new LinkedHashSet<>();
			final SAMFileHeader header = SamReaderFactory.
				makeDefault().
				validationStringency(ValidationStringency.SILENT).
				getFileHeader(path);
			for(final SAMReadGroupRecord rg:header.getReadGroups())
				{
				final StringBuilder line = new StringBuilder();
				boolean found_one = false;
				
				for(final String rgAttribute: attributes) {
					final String sn = rg.getAttribute( rgAttribute );
					if (sn != null) found_one = true;
					line.append(sn==null?"":sn).append("\t");
					}
				if(!found_one) continue;
				line.append(pathStr);
				samples.add(line.toString());
				}
			for(final String sampleLine: samples) lines.add(sampleLine);
			}
		else  if(pathStr.endsWith(FileExtensions.VCF) ||
			pathStr.endsWith(FileExtensions.COMPRESSED_VCF) ||
			pathStr.endsWith(FileExtensions.BCF)
			)
			{
			VCFFileReader vcfFileReader = null;
			try {
				vcfFileReader = new VCFFileReader(path, false);
				final VCFHeader header = vcfFileReader.getFileHeader()
				for(final String sn: header.getSampleNamesInOrder())
					{
					lines.add(sn + "\t" + pathStr);
					}
				}
			finally {
				if (vcfFileReader!=null) vcfFileReader.close(); 
				}
			}
		else
			{
			throw new IllegalArgumentException("file format is not supported. \"" + path + "\" is not a VCF/BAM/SAM/CRAM file.");
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
    SamplesSplitter options(Map opts) {
        super.options(opts)
        if(opts.attribute) this.rgAttributes = opts.attribute;
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
