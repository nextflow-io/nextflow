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

package nextflow.htsjdk

import java.nio.charset.Charset
import java.nio.file.Paths
import java.nio.file.Files
import java.nio.file.Path
import java.io.File

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j


import htsjdk.variant.utils.SAMSequenceDictionaryExtractor
import htsjdk.samtools.SAMSequenceDictionary
import htsjdk.samtools.SAMSequenceRecord
import htsjdk.samtools.util.FileExtensions
import htsjdk.samtools.SamReaderFactory
import htsjdk.samtools.SAMFileHeader
import htsjdk.samtools.SAMReadGroupRecord
import htsjdk.samtools.ValidationStringency
import htsjdk.variant.vcf.VCFFileReader
import htsjdk.variant.vcf.VCFHeader

/**
 * Utilities for htsjdk
 *
 * @author Pierre Lindenbaum PhD Institut du Thorax Nantes France.
 */
@Slf4j
@CompileStatic
class HtsjdkUtils  {
     /* base class working on BAM/VCF/etc.. files */
    private static abstract class AbstractPathClosure
	{
	protected abstract List applyPath(Path p);

	public List apply(obj) {
		if( obj == null )
			{
			throw new IllegalArgumentException("Cannot apply on null Path.");
			}
		else if( obj instanceof Path) {
			return applyPath((Path)obj);
			}
		else if( obj instanceof String) {
			return applyPath(Paths.get((String)obj) );
			}
		else if( obj instanceof File) {
			return applyPath(((File)obj).toPath());
			}
		else
			{
			throw new IllegalArgumentException("Cannot create Path from "+obj.getClass());
			}
		}
	}

   /** extracts a SequenceDictionary from a path */
   private static class DictionaryClosure extends AbstractPathClosure
	{
	@Override
	protected List applyPath(final Path path) {
		final SAMSequenceDictionary dict = SAMSequenceDictionaryExtractor.extractDictionary(path);
		if(dict==null || dict.isEmpty()) return Collections.emptyList();
		def rows = new ArrayList<>(dict.size());
		for(final SAMSequenceRecord ssr : dict.getSequences() ) {
			def attributes = new HashMap(5);
			attributes.put(	"index" , ssr.getSequenceIndex());
			attributes.put( SAMSequenceRecord.SEQUENCE_NAME_TAG , ssr.getSequenceName());
			attributes.put( SAMSequenceRecord.SEQUENCE_LENGTH_TAG , ssr.getSequenceLength());
			if(ssr.getAssembly()!=null) attributes.put( SAMSequenceRecord.ASSEMBLY_TAG , ssr.getAssembly());
			if(ssr.getSpecies()!=null) attributes.put( SAMSequenceRecord.SPECIES_TAG , ssr.getSpecies());
			def row = [attributes,path]
			rows.add(row);
			}
		return rows;
		}
	}
   
    /** extracts the samples (VCF) or the read groups (BAM/SAM) */
   private static class SampleClosure extends AbstractPathClosure
	{
	private List fromSamOrBam(final Path path) {
		def rows = new ArrayList<>();
		final SAMFileHeader header = SamReaderFactory.
				makeDefault().
				validationStringency(ValidationStringency.SILENT).
				getFileHeader(path);
		for(final SAMReadGroupRecord rg:header.getReadGroups())
			{
			rows.add([ rg.getAttributes() , path]);
			}
		return rows;
		}

	private List fromVcf(final Path path) {
		VCFFileReader vcfFileReader = null;
		def rows = new ArrayList<>();
		try {
			vcfFileReader = new VCFFileReader(path, false);
			final VCFHeader header = vcfFileReader.getFileHeader()
			for(final String sn: header.getSampleNamesInOrder())
				{
				def attributes = new HashMap(1);
				attributes.put( SAMReadGroupRecord.READ_GROUP_SAMPLE_TAG, sn);
				rows.add([ attributes , path]);
				}
			}
		finally {
			if (vcfFileReader!=null) vcfFileReader.close(); 
			}
		return rows;
		}

	@Override
	protected List applyPath(final Path path) {
		final String pathStr = path.toString();
		if( pathStr.endsWith(FileExtensions.SAM) ||
			pathStr.endsWith(FileExtensions.BAM) 
			)
			{
			return fromSamOrBam(path);
			}
		else  if(pathStr.endsWith(FileExtensions.VCF) ||
			pathStr.endsWith(FileExtensions.COMPRESSED_VCF) ||
			pathStr.endsWith(FileExtensions.BCF)
			)
			{
			return fromVcf(path);
			}
		else
			{
			throw new IllegalArgumentException("Cannot extract samples from "+path);
			}
		}
	}

   public static Closure dictionary() {
	final DictionaryClosure extractor = new DictionaryClosure();
	return { obj-> extractor.apply(obj) };
	}
   public static Closure samples() {
	final SampleClosure extractor = new SampleClosure();
	return { obj-> extractor.apply(obj) };
	}
}
