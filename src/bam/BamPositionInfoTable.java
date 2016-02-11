package bam;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.MappedFragment;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotation.predicate.SecondReadFilter;
import guttmanlab.core.annotationcollection.AbstractAnnotationCollection;
import guttmanlab.core.annotationcollection.AnnotationCollection;
import guttmanlab.core.annotationcollection.BAMFragmentCollectionFactory;
import guttmanlab.core.annotationcollection.BAMPairedFragmentCollection;
import guttmanlab.core.annotationcollection.FilteredIterator;
import guttmanlab.core.util.CommandLineParser;
import guttmanlab.core.util.CountLogger;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Map;
import java.util.TreeMap;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import net.sf.samtools.util.CloseableIterator;

public class BamPositionInfoTable {
	
	private static Logger logger = Logger.getLogger(BamPositionInfoTable.class.getName());
	
	/**
	 * Get counts of fragments starting at each position in the region
	 * Includes positions outside the regions where overlapping fragments start
	 * @param data Alignment data
	 * @param region Region consisting of coordinates to get
	 * @param enforceSameOrientation Only count fragments with same orientation as region
	 * @return Map of position to number of fragments whose 5' end is at the position
	 */
	public static Map<Integer, Integer> getStartPosCounts(AbstractAnnotationCollection<? extends MappedFragment> data, 
			Annotation region, boolean enforceSameOrientation) {
		logger.debug("");
		logger.debug("");
		logger.debug("Getting start position counts for " + region.toBED());
		CloseableIterator<? extends MappedFragment> iter = data.sortedIterator(region, false);
		Map<Integer, Integer> rtrn = new TreeMap<Integer, Integer>();
		Strand regionOrientation = region.getOrientation();
		if(enforceSameOrientation) {
			if(!regionOrientation.equals(Strand.POSITIVE) && !regionOrientation.equals(Strand.NEGATIVE)) {
				throw new IllegalArgumentException("Can't enforce strand unless region strand is positive or negative");
			}
		}
		int numSkipped = 0;
		while(iter.hasNext()) {
			MappedFragment fragment = iter.next();
			logger.debug("Next fragment\t" + fragment.toBED());
			Strand orientation = fragment.getOrientation();
			if(enforceSameOrientation && !orientation.equals(regionOrientation)) {
				numSkipped++;
				continue;
			}
			int fragStart = -1;
			switch(orientation) {
			case BOTH:
				iter.close();
				throw new IllegalArgumentException("Strand must be positive or negative");
			case INVALID:
				iter.close();
				throw new IllegalArgumentException("Strand must be positive or negative");
			case NEGATIVE:
				fragStart = fragment.getReferenceEndPosition() - 1;
				break;
			case POSITIVE:
				fragStart = fragment.getReferenceStartPosition();
				break;
			case UNKNOWN:
				iter.close();
				throw new IllegalArgumentException("Strand must be positive or negative");
			default:
				iter.close();
				throw new IllegalArgumentException("Strand must be positive or negative");
			}
			
			Integer key = Integer.valueOf(fragStart);
			if(!rtrn.containsKey(key)) {
				rtrn.put(key, Integer.valueOf(0));
			}
			rtrn.put(key, Integer.valueOf(rtrn.get(key).intValue() + 1));
			logger.debug("Start pos: " + fragStart + " key,value: " + key.toString() + "," + rtrn.get(key).toString());
		}
		iter.close();
		if(enforceSameOrientation) logger.info("Skipped " + numSkipped + " records that did not match orientation for " + region.getName());
		return rtrn;
	}
	
	/**
	 * Append counts to table
	 * @param writer Existing writer
	 * @param chr Chromosome name
	 * @param counts Map of position to count
	 * @throws IOException
	 */
	private static void append(BufferedWriter writer, String chr, Map<Integer, Integer> counts) throws IOException {
		for(Integer pos : counts.keySet()) {
			writer.write(chr + "\t" + pos.toString() + "\t" + counts.get(pos).toString() + "\n");
		}
	}
	
	private static void updateCounts(Map<String, Map<Integer, Integer>> countsByChr, String chr, Map<Integer, Integer> counts) {
		if(!countsByChr.containsKey(chr)) {
			countsByChr.put(chr, new TreeMap<Integer, Integer>());
		}
		for(Integer pos : counts.keySet()) {
			if(!countsByChr.get(chr).containsKey(pos)) {
				countsByChr.get(chr).put(pos, Integer.valueOf(0));
			}
			int oldVal = countsByChr.get(chr).get(pos).intValue();
			int toAdd = counts.get(pos).intValue();
			countsByChr.get(chr).put(pos, Integer.valueOf(oldVal + toAdd));
		}
	}
	
	private static void writeCountTable(String bamFile, String regionBed, String chrSizeFile, String outTable, 
			boolean enforceSameOrientation, boolean forceSingleEnd, boolean secondReadOnly) throws IOException {
		AbstractAnnotationCollection<? extends MappedFragment> data = BAMFragmentCollectionFactory.createFromBam(bamFile, forceSingleEnd);
		if(secondReadOnly) data.addFilter(new SecondReadFilter());
		BufferedWriter writer = new BufferedWriter(new FileWriter(outTable));
		AnnotationCollection<Gene> regions = BEDFileIO.loadFromFile(regionBed, chrSizeFile);
		int numGenes = regions.getNumAnnotations();
		CountLogger cl = new CountLogger(numGenes, 100);
		CloseableIterator<Gene> iter = regions.sortedIterator();
		Map<String, Map<Integer, Integer>> countsByChr = new TreeMap<String, Map<Integer, Integer>>();
		while(iter.hasNext()) {
			Gene region = iter.next();
			cl.advance();
			String chr = region.getReferenceName();
			Map<Integer, Integer> counts = getStartPosCounts(data, region, enforceSameOrientation);
			updateCounts(countsByChr, chr, counts);
		}
		for(String chr : countsByChr.keySet()) {
			append(writer, chr, countsByChr.get(chr));
		}
		iter.close();
		writer.close();
	}
	
	
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-bam", "Bam file", true);
		p.addStringArg("-bed", "Bed file of regions. Will write count for every position in blocks", true);
		p.addStringArg("-c", "Chromosome size file", true);
		p.addStringArg("-o", "Output table", true);
		p.addBooleanArg("-eo", "Only count fragments mapped in same orientation as gene they overlap", true);
		p.addBooleanArg("-se", "Force single end interpretation of mappings", true);
		p.addBooleanArg("-d", "Debug logging", false, false);
		p.addBooleanArg("-r2", "Count read 2 only", true);
		p.parse(args);
		if(p.getBooleanArg("-d")) {
			BAMPairedFragmentCollection.logger.setLevel(Level.DEBUG);
			FilteredIterator.logger.setLevel(Level.DEBUG);
			logger.setLevel(Level.DEBUG);
		}
		String bamFile = p.getStringArg("-bam");
		String regionBed = p.getStringArg("-bed");
		String chrSizeFile = p.getStringArg("-c");
		String outTable = p.getStringArg("-o");
		boolean enforceSameOrientation = p.getBooleanArg("-eo");
		boolean forceSingleEnd = p.getBooleanArg("-se");
		boolean read2only = p.getBooleanArg("-r2");
		
		writeCountTable(bamFile, regionBed, chrSizeFile, outTable, enforceSameOrientation, forceSingleEnd, read2only);
		
		logger.info("");
		logger.info("All done.");
		
	}
	
}
