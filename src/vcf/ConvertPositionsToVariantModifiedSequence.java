package vcf;

import guttmanlab.core.util.CommandLineParser;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

/**
 * Convert positions in a VCF file to positions on a modified reference sequence that incorporates all variants in the sequence
 * @author prussell
 *
 */
public final class ConvertPositionsToVariantModifiedSequence {
	
	// Prohibit instantiation
	private ConvertPositionsToVariantModifiedSequence() {}

	/**
	 * Convert positions to positions on a modified reference sequence that incorporates all variants in the input set
	 * These should be all the variants for their reference sequence
	 * Throws exception if any records come from different reference sequences
	 * @param origRecords All the records for some reference sequence
	 * @return Records converted to positions on a modified sequence where all the input variants have been incorporated
	 */
	@SuppressWarnings("unused")
	private static TreeSet<VCFRecord> convertPositions(TreeSet<VCFRecord> origRecords) {
		return convertPositions(origRecords, "");
	}
	
	/**
	 * Convert a whole VCF file and write new file
	 * @param origVcf Original VCF file
	 * @param outVcf Output VCF file
	 * @param appendToRefName String to append to reference sequence name for each returned record
	 * @throws IOException
	 */
	private static void convertFile(String origVcf, String outVcf, String appendToRefName) throws IOException {
		Map<String, TreeSet<VCFRecord>> recordsByChr = fromVcfFile(origVcf);
		FileWriter w = new FileWriter(outVcf);
		BufferedReader b = new BufferedReader(new FileReader(origVcf));
		while(b.ready()) { // Copy header section
			String line = b.readLine();
			if(line.startsWith("#")) {
				w.write(line + "\n");
			} else break;
		}
		for(String chr : recordsByChr.keySet()) {
			TreeSet<VCFRecord> converted = convertPositions(recordsByChr.get(chr), appendToRefName);
			for(VCFRecord record : converted) {
				w.write(record.toString() + "\n");
			}
		}
		b.close();
		w.close();
	}
	
	/**
	 * Read records from VCF file and organize by reference name
	 * @param file VCF file
	 * @return Map of reference sequence name to sorted set of all records on the reference sequence
	 * @throws IOException
	 */
	private static Map<String, TreeSet<VCFRecord>> fromVcfFile(String file) throws IOException {
		BufferedReader b = new BufferedReader(new FileReader(file));
		Map<String, TreeSet<VCFRecord>> rtrn = new TreeMap<String, TreeSet<VCFRecord>>();
		String header = null;
		while(b.ready()) {
			String line = b.readLine();
			if(line.startsWith("##")) continue;
			if(line.startsWith("#")) {
				header = line;
				break;
			}
		}
		if(!b.ready()) {
			b.close();
			throw new IllegalArgumentException("Nothing after header line");
		}
		while(b.ready()) {
			String line = b.readLine();
			if(line.substring(0, 1).equals("#")) {
				b.close();
				throw new IllegalArgumentException("Comment line after header line");
			}
			VCFRecord record = new VCFRecord(line, header);
			String chr = record.getChrom();
			if(!rtrn.containsKey(chr)) {
				rtrn.put(chr, new TreeSet<VCFRecord>());
			}
			rtrn.get(chr).add(record);
		}
		b.close();
		return rtrn;
	}
		
	/**
	 * Convert positions to positions on a modified reference sequence that incorporates all variants in the input set
	 * These should be all the variants for their reference sequence
	 * Throws exception if any records come from different reference sequences
	 * @param origRecords All the records for some reference sequence
	 * @param appendToRefName String to append to reference sequence name for each returned record
	 * @return Records converted to positions on a modified sequence where all the input variants have been incorporated
	 */
	private static TreeSet<VCFRecord> convertPositions(TreeSet<VCFRecord> origRecords, String appendToRefName) {
		Iterator<VCFRecord> descIter = origRecords.descendingIterator();
		if(origRecords.isEmpty()) {
			throw new IllegalArgumentException("Set of records is empty");
		}
		VCFRecord last = descIter.next();
		String chr = last.getChrom();
		String modifiedChr = chr + appendToRefName;

		TreeSet<VCFRecord> inProgress = new TreeSet<VCFRecord>();
		inProgress.add(VCFRecord.modifiedCopy(last, modifiedChr, last.getOneBasedPos()));
		while(descIter.hasNext()) {
			VCFRecord record = descIter.next();
			Collection<VCFRecord> modifiedForThisRecord = new TreeSet<VCFRecord>();
			modifiedForThisRecord.add(VCFRecord.modifiedCopy(record, modifiedChr, record.getOneBasedPos()));
			int offset = record.getAlternateAllele().getLength() - record.getRefAllele().getLength();
			for(VCFRecord prev : inProgress) {
				modifiedForThisRecord.add(VCFRecord.modifiedCopy(prev, modifiedChr, prev.getOneBasedPos() + offset));
			}
			inProgress.clear();
			inProgress.addAll(modifiedForThisRecord);
		}
		return inProgress;
	}
	
	
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-i", "Input VCF", true);
		p.addStringArg("-o", "Output VCF", true);
		p.addStringArg("-a", "Append string to reference names", false, "");
		p.parse(args);
		convertFile(p.getStringArg("-i"), p.getStringArg("-o"), p.getStringArg("-a"));
		
	}
	
	
}
