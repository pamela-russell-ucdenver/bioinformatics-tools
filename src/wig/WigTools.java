/**
 * 
 */
package wig;

import guttmanlab.core.util.CommandLineParser;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.log4j.Logger;

import broad.core.sequence.Sequence;
import broad.core.sequence.TranscribedRegions;
import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Gene;
import nextgen.core.readers.WigReader;
import nextgen.core.utils.CountLogger;
import nextgen.core.writers.WigWriter;

/**
 * @author prussell
 *
 */
public class WigTools {

	private WigReader wigReader;
	private Map<String, TreeMap<Integer,Double>> wigData;
	private static Logger logger = Logger.getLogger(WigTools.class.getName());
	private Map<String, Double> baseCounts;
	private Map<String, Double> dinucleotideCounts;
	private TranscribedRegions transcribedRegions;
	private static String[] bases = {"A", "C", "G", "T", "N"};
	
	/**
	 * @param genomeFasta Genome fasta file
	 * @param bedFile Bed file of genes to get orientations from
	 * @throws IOException
	 */
	public WigTools(String genomeFasta, String bedFile) throws IOException {
		this(null, genomeFasta, bedFile);
		logger.warn("No input wig file provided");
	}
	
	/**
	 * @param wigFile Input wig file
	 * @param genomeFasta Genome fasta file
	 * @param bedFile Bed file of genes to get orientations from
	 * @throws IOException
	 */
	public WigTools(String wigFile, String genomeFasta, String bedFile) throws IOException {
		transcribedRegions = new TranscribedRegions(genomeFasta, bedFile);
		// Read wig file
		if(wigFile != null) {
			wigReader = new WigReader(wigFile);
			wigData = wigReader.getAllValues();
		}
		baseCounts = new TreeMap<String, Double>();
		dinucleotideCounts = new TreeMap<String, Double>();
	}
	
	/**
	 * Shift all values in wig file to another position
	 * Shift is specified by an offset, in direction of transcription
	 * The new position is obtained by calling TranscribedRegions.shiftPosition()
	 * The value of the original position is assigned to the shifted position in the return structure
	 * @param offset Offset along direction of transcription (can be negative)
	 * @return Shifted values by chromosome
	 */
	private Map<String, TreeMap<Integer,Double>> shiftWig(int offset) {
		logger.info("Getting track with all values shifted along direction of transcription...");
		Map<String, TreeMap<Integer, Double>> rtrn = new TreeMap<String, TreeMap<Integer, Double>>();
		for(String chr : wigData.keySet()) {
			logger.info(chr);
			TreeMap<Integer, Double> shiftedCountsThisChr = new TreeMap<Integer, Double>();
			for(Integer origPos : wigData.get(chr).keySet()) {
				try {
					Integer shiftedPos = Integer.valueOf(transcribedRegions.shiftPosition(chr, origPos.intValue(), offset));
					Double shiftedVal = wigData.get(chr).get(origPos);
					shiftedCountsThisChr.put(shiftedPos, shiftedVal);
				} catch (IllegalArgumentException e) {
					logger.debug("Skipping original position " + chr + " " + origPos.toString() + ": " + e.getMessage());
				} catch (UnsupportedOperationException e) {
					logger.debug("Skipping original position " + chr + " " + origPos.toString() + ": " + e.getMessage());
				}
			}
			rtrn.put(chr, shiftedCountsThisChr);
		}
		logger.info("Done getting shifted track.");
		return rtrn;
	}
	
	/**
	 * Write shifted values returned by shiftWig()
	 * @param outFile Output wig file
	 * @param offset Offset along direction of transcription (can be negative)
	 * @throws IOException
	 */
	private void writeShiftedWig(String outFile, int offset) throws IOException {
		Map<String, TreeMap<Integer, Double>> shiftedWig = shiftWig(offset);
		FileWriter w = new FileWriter(outFile);
		logger.info("Writing shifted wig file to " + outFile + "...");
		for(String chr : shiftedWig.keySet()) {
			logger.info(chr);
			WigWriter.write(w, chr, shiftedWig.get(chr), false);
		}
		w.close();
		logger.info("Done writing shifted file.");
	}
	
	/**
	 * Get a track for each possible transcribed nucleotide
	 * The underlying wig data is divided into tracks for each nucleotide
	 * Direction of transcription is discovered based on the annotation provided
	 * If position is not transcribed or direction of transcription is ambiguous, position is not included
	 * @return Map of nucleotide to chromosome-based count data
	 * @throws IOException 
	 */
	private Map<String, Map<String, TreeMap<Integer, Double>>> getTracksByNucleotide() throws IOException {
		logger.info("Getting separate tracks for each nucleotide...");
		// Initialize data structure
		Map<String, Map<String, TreeMap<Integer, Double>>> rtrn = new TreeMap<String, Map<String, TreeMap<Integer, Double>>>();
		for(int i=0; i<bases.length; i++) {
			Map<String, TreeMap<Integer, Double>> trackThisBase = new TreeMap<String, TreeMap<Integer, Double>>();
			for(String chr : wigData.keySet()) {
				TreeMap<Integer, Double> trackThisChr = new TreeMap<Integer, Double>();
				trackThisBase.put(chr, trackThisChr);
			}
			rtrn.put(bases[i], trackThisBase);
		}
		// Make tracks
		for(String chr : wigData.keySet()) {
			logger.info(chr);
			for(Integer pos : wigData.get(chr).keySet()) {
				Double val = wigData.get(chr).get(pos);
				char[] base = {transcribedRegions.getTranscribedBase(chr, pos.intValue())};
				String baseString = new String(base);
				rtrn.get(baseString).get(chr).put(pos, val);
			}
		}
		logger.info("Done getting separate tracks for each nucleotide.");
		return rtrn;
	}
	
	/**
	 * Get separate tracks for each nucleotide and write to separate wig files
	 * @param outFilePrefix Prefix for output wig files
	 * @throws IOException
	 */
	private void writeTracksByNuceotide(String outFilePrefix) throws IOException{
		logger.info("Writing separate nucleotide tracks...");
		Map<String, Map<String, TreeMap<Integer, Double>>> tracksByNuc = getTracksByNucleotide();
		for(int i=0; i<bases.length; i++) {
			String base = bases[i];
			String outFile = outFilePrefix + "_" + base + ".wig";
			FileWriter w = new FileWriter(outFile);
			for(String chr : tracksByNuc.get(base).keySet()) {
				WigWriter.write(w, chr, tracksByNuc.get(base).get(chr), false);
			}
			w.close();
		}
		logger.info("Done writing separate tracks.");
	}
	
	
	/**
	 * Add values from two wig files and write merged wig file
	 * @param wig1 Wig file 1
	 * @param wig2 Wig file 2
	 * @param outFile Output file to write
	 * @throws IOException
	 */
	private static void mergeWigFiles(String wig1, String wig2, String outFile) throws IOException {
		logger.info("Merging wig files " + wig1 + " and " + wig2 + " and writing to file " + outFile + "...");
		WigReader r1 = new WigReader(wig1);
		WigReader r2 = new WigReader(wig2);
		Map<String, TreeMap<Integer, Double>> vals1 = r1.getAllValues();
		Map<String, TreeMap<Integer, Double>> vals2 = r2.getAllValues();
		Collection<String> allChrs = new TreeSet<String>();
		allChrs.addAll(vals1.keySet());
		allChrs.addAll(vals2.keySet());
		FileWriter w = new FileWriter(outFile);
		for(String chr : allChrs) {
			logger.info(chr);
			if(!vals1.containsKey(chr)) {
				WigWriter.write(w, chr, vals2.get(chr), false);
				continue;
			}
			if(!vals2.containsKey(chr)) {
				WigWriter.write(w, chr, vals1.get(chr), false);
				continue;
			}
			TreeMap<Integer, Double> valsThisChr = new TreeMap<Integer, Double>();
			Collection<Integer> positions = new TreeSet<Integer>();
			positions.addAll(vals1.get(chr).keySet());
			positions.addAll(vals2.get(chr).keySet());
			for(Integer pos : positions) {
				double val = 0;
				if(vals1.get(chr).containsKey(pos)) {
					val += vals1.get(chr).get(pos).doubleValue();
				}
				if(vals2.get(chr).containsKey(pos)) {
					val += vals2.get(chr).get(pos).doubleValue();
				}
				valsThisChr.put(pos, Double.valueOf(val));
			}
			WigWriter.write(w, chr, valsThisChr, false);
		}
		w.close();
		logger.info("Done writing merged wig file.");
	}
	
	
	private static void writeCountsToTable(Map<String, Double> counts, String outFile) throws IOException {
		FileWriter w = new FileWriter(outFile);
		double total = 0;
		for(String base : counts.keySet()) {
			total += counts.get(base).doubleValue();
		}
		for(String base : counts.keySet()) {
			double count = counts.get(base).doubleValue();
			w.write(base + "\t" + count + "\t" + count / total + "\n");
		}
		w.close();		
	}
	
	/**
	 * Write base counts to table
	 * @param outFile Output file for table
	 * @throws IOException
	 */
	private void writeBaseCountsToTable(String outFile) throws IOException {
		if(baseCounts.isEmpty()) {
			computeBaseCounts();
		}
		writeCountsToTable(baseCounts, outFile);
	}
	
	/**
	 * Write dinucleotide counts to table
	 * @param firstPosRelative First position of dinucleotide relative to the wig position, in transcriptome coordinates
	 * @param outFile Output file for table
	 * @throws IOException
	 */
	private void writeDinucCountsToTable(int firstPosRelative, String outFile) throws IOException {
		if(dinucleotideCounts.isEmpty()) {
			computeDinucleotideCounts(firstPosRelative);
		}
		writeCountsToTable(dinucleotideCounts, outFile);
	}
	
	
	
	/**
	 * Compute counts for each base
	 * For each position in wig file, first identify the transcribed base at that position using gene annotation
	 * Next increment count for that base by adding the value from the wig file at the position
	 * @throws IOException 
	 */
	private void computeBaseCounts() throws IOException {
		logger.info("Computing base counts...");
		double[] counts = new double[10];
		for(int i = 0; i < counts.length; i++) {
			counts[i] = 0;
		}
		for(String chr : wigData.keySet()) {
			logger.info(chr);
			for(Integer pos : wigData.get(chr).keySet()) {
				double value = wigData.get(chr).get(pos).doubleValue();
				char base = transcribedRegions.getTranscribedBase(chr, pos.intValue());
				switch(base) {
					case 'A':  counts[Sequence.SHORT_ENCODED_A] += value;
					break;
					case 'C':  counts[Sequence.SHORT_ENCODED_C] += value;
					break;
					case 'G':  counts[Sequence.SHORT_ENCODED_G] += value;
					break;
					case 'T':  counts[Sequence.SHORT_ENCODED_T] += value;
					break;
					case 'a':  counts[Sequence.SHORT_ENCODED_a] += value;
					break;
					case 'c':  counts[Sequence.SHORT_ENCODED_c] += value;
					break;
					case 'g':  counts[Sequence.SHORT_ENCODED_g] += value;
					break;
					case 't':  counts[Sequence.SHORT_ENCODED_t] += value;
					break;
					case 'N':  counts[Sequence.SHORT_ENCODED_N] += value;
					break;
					default: throw new IllegalArgumentException("Base " + base + " not valid.");
				}
			}
		}
		baseCounts.put("A", Double.valueOf(counts[Sequence.SHORT_ENCODED_A]));
		baseCounts.put("C", Double.valueOf(counts[Sequence.SHORT_ENCODED_C]));
		baseCounts.put("G", Double.valueOf(counts[Sequence.SHORT_ENCODED_G]));
		baseCounts.put("T", Double.valueOf(counts[Sequence.SHORT_ENCODED_T]));
		baseCounts.put("a", Double.valueOf(counts[Sequence.SHORT_ENCODED_a]));
		baseCounts.put("c", Double.valueOf(counts[Sequence.SHORT_ENCODED_c]));
		baseCounts.put("g", Double.valueOf(counts[Sequence.SHORT_ENCODED_g]));
		baseCounts.put("t", Double.valueOf(counts[Sequence.SHORT_ENCODED_t]));
		baseCounts.put("N", Double.valueOf(counts[Sequence.SHORT_ENCODED_N]));
		logger.info("Done computing base counts.");
	}
	
	/**
	 * Compute counts for each dinucleotide
	 * Converts all bases to upper case
	 * @param firstPosRelative First position of dinucleotide relative to the wig position, in transcriptome coordinates
	 * @throws IOException
	 */
	private void computeDinucleotideCounts(int firstPosRelative) throws IOException {
		
		logger.info("Computing dinucleotide counts...");
		dinucleotideCounts.clear();
		
		for(String chr : wigData.keySet()) {
			logger.info(chr);
			for(Integer pos : wigData.get(chr).keySet()) {
				double value = wigData.get(chr).get(pos).doubleValue();
				try {
					char firstBase = transcribedRegions.getTranscribedBase(chr, transcribedRegions.shiftPosition(chr, pos.intValue(), firstPosRelative));
					char secondBase = transcribedRegions.getTranscribedBase(chr, transcribedRegions.shiftPosition(chr, pos.intValue(), firstPosRelative + 1));
					char[] dinuc = new char[2];
					dinuc[0] = firstBase;
					dinuc[1] = secondBase;
					String dinucString = new String(dinuc).toUpperCase();
					if(!dinucleotideCounts.containsKey(dinucString)) {
						dinucleotideCounts.put(dinucString, Double.valueOf(value));
						continue;
					}
					dinucleotideCounts.put(dinucString, Double.valueOf(dinucleotideCounts.get(dinucString).doubleValue() + value));
				} catch(UnsupportedOperationException e) {
					logger.debug("Skipping " + chr + ":" + pos.toString());
					continue;
				} catch(IllegalArgumentException e) {
					logger.debug("Skipping " + chr + ":" + pos.toString());
					continue;
				}
			}			
		}
		logger.info("Done computing dinucleotide counts.");
	}
	
	/**
	 * Write tracks of transcribed nucleotides
	 * One track for each of A, C, G, T
	 * Wig value is 1 if that base is transcribed at the given position
	 * @param outWigFilePrefix Prefix for output wig tracks
	 * @throws IOException
	 */
	private void writeBinaryNucleotideTracks(String outWigFilePrefix) throws IOException {
		String fileA = outWigFilePrefix + "_A.wig";
		String fileC = outWigFilePrefix + "_C.wig";
		String fileG = outWigFilePrefix + "_G.wig";
		String fileT = outWigFilePrefix + "_T.wig";
		logger.info("Writing binary transcribed nucleotide tracks to files " + fileA + ", " + fileC + ", " + fileG + ", " + fileT + "...");
		FileWriter writerA = new FileWriter(fileA);
		FileWriter writerC = new FileWriter(fileC);
		FileWriter writerG = new FileWriter(fileG);
		FileWriter writerT = new FileWriter(fileT);
		Map<String, Collection<Gene>> genes = transcribedRegions.getGenes();
		for(String chr : genes.keySet()) {
			logger.info(chr);
			// Initialize the four tracks
			TreeMap<Integer, Double> transcribedA = new TreeMap<Integer, Double>();
			TreeMap<Integer, Double> transcribedC = new TreeMap<Integer, Double>();
			TreeMap<Integer, Double> transcribedG = new TreeMap<Integer, Double>();
			TreeMap<Integer, Double> transcribedT = new TreeMap<Integer, Double>();
			// First list all postions contained in any exon
			Collection<Integer> allPositions = new TreeSet<Integer>();
			for(Gene gene : genes.get(chr)) {
				for(Annotation exon : gene.getBlocks()) {
					for(int i = exon.getStart(); i < exon.getEnd(); i++) {
						allPositions.add(Integer.valueOf(i));
					}
				}
			}
			// Now get transcribed nucleotide for each position and add to tracks
			int numPositions = allPositions.size();
			CountLogger c = new CountLogger(numPositions,10);
			for(Integer pos : allPositions) {
				c.advance();
				char base = transcribedRegions.getTranscribedBase(chr, pos.intValue());
				switch(base) {
				case 'A':
					transcribedA.put(pos, Double.valueOf(1));
					break;
				case 'a':
					transcribedA.put(pos, Double.valueOf(1));
					break;
				case 'C':
					transcribedC.put(pos, Double.valueOf(1));
					break;
				case 'c':
					transcribedC.put(pos, Double.valueOf(1));
					break;
				case 'G':
					transcribedG.put(pos, Double.valueOf(1));
					break;
				case 'g':
					transcribedG.put(pos, Double.valueOf(1));
					break;
				case 'T':
					transcribedT.put(pos, Double.valueOf(1));
					break;
				case 't':
					transcribedT.put(pos, Double.valueOf(1));
					break;
				default:
					break;
				}
			}
			// Write tracks
			WigWriter.write(writerA, chr, transcribedA, false);
			WigWriter.write(writerC, chr, transcribedC, false);
			WigWriter.write(writerG, chr, transcribedG, false);
			WigWriter.write(writerT, chr, transcribedT, false);
		}
		writerA.close();
		writerC.close();
		writerG.close();
		writerT.close();
	}
	

	/**
	 * Intersect a gene with the wig data
	 * @param gene The gene
	 * @param minWigValue Minimum wig value to include in intersection
	 * @param identifier Identifier to append to name of returned gene
	 * @return A gene whose blocks consist of contiguous positions with wig value over the min, intersected with the gene
	 */
	private Gene getWigOverlappingPositionsAsGene(Gene gene, double minWigValue, String identifier) {
		String chr = gene.getChr();
		// Get all wig values for positions overlapping exons of the gene
		TreeMap<Integer, Double> exonicWigValues = new TreeMap<Integer, Double>();
		for(Annotation exon : gene.getExonSet()) {
			exonicWigValues.putAll(wigReader.getValues(chr, exon.getStart(), exon.getEnd()));
		}
		if(exonicWigValues.isEmpty()) {
			logger.debug("Gene " + gene.getName() + " has no overlapping wig positions with value >" + minWigValue + ". Returning null.");
			return null;
		}
		// Make contiguous positions with wig value over min into blocks
		Iterator<Integer> posIter = exonicWigValues.keySet().iterator();
		Collection<Gene> contiguousBlocks = new TreeSet<Gene>();
		int start = posIter.next().intValue();
		int end = start;
		while(posIter.hasNext()) {
			Integer n = posIter.next();
			if(wigData.get(chr).get(n).doubleValue() < minWigValue) {
				continue;
			}
			int next = n.intValue();
			if(next == end + 1)	{
				end = next;
			} else {
				Gene newBlock = new Gene(chr, start, end + 1);
				contiguousBlocks.add(newBlock);
				start = next;
				end = next;
			}
		}
		Gene newBlock = new Gene(chr, start, end + 1);
		contiguousBlocks.add(newBlock);
		// Overlap the blocks with the gene
		Gene rtrn = gene.getOverlap(contiguousBlocks);
		rtrn.setName(gene.getName() + "_" + identifier);
		return rtrn;
	}
	
	/**
	 * Intersect all genes with the wig data and write to a bed file
	 * The bed lines are genes whose blocks consist of contiguous positions with wig value over the min, intersected with each gene
	 * @param outBed Output bed file
	 * @param minWigValue Minimum wig value to include in intersection
	 * @param identifier Identifier to append to names of genes written to bed file
	 * @throws IOException
	 */
	private void writeWigOverlappingPositionsAsBed(String outBed, double minWigValue, String identifier) throws IOException {
		FileWriter w = new FileWriter(outBed);
		logger.info("Writing genes intersected with wig data to file " + outBed + "...");
		for(String chr : transcribedRegions.getGenes().keySet()) {
			logger.info(chr);
			int numGenes = transcribedRegions.getGenes().get(chr).size();
			CountLogger c = new CountLogger(numGenes, 10);
			for(Gene gene : transcribedRegions.getGenes().get(chr)) {
				c.advance();
				Gene intersect = getWigOverlappingPositionsAsGene(gene, minWigValue, identifier);
				if(intersect == null) continue;
				w.write(intersect.toBED() + "\n");
			}
		}
		w.close();
		logger.info("Done writing file.");
	}
	
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		String description = "\n***** WIG TOOLS *****\n";
		description += "Task 1: write table of total base counts (requires -w, -g, -b, -ot)\n";
		description += "Task 2: write tracks for each nucleotide (requires -w, -g, -b, -on)\n";
		description += "Task 3: write shifted wig file (requires -w, -g, -b, -of, -os)\n";
		description += "Task 4: merge two wig files (requires -w, -w2, -om)\n";
		description += "Task 5: write binary wig tracks of transcribed nucleotides at all exonic positions (requires -g, -b, -ob)\n";
		description += "Task 6: write bed file of genes whose blocks consist of contiguous positions with wig value over a min value, intersected with each gene (requires -w, -g, -b, -mw, -id, -oi)\n";
		description += "Task 7: write table of total dinucleotide counts (requires -w, -g, -b, -otd)\n";
		
		p.setProgramDescription(description);
		p.addStringArg("-g", "For task 1, 2, 3, 5 or 7: genome fasta", false, null);
		p.addStringArg("-b", "For task 1, 2, 3, 5 or 7: bed annotation", false, null);
		p.addStringArg("-w2", "For task 4: other wig file to merge", false, null);
		p.addStringArg("-w", "For task 1, 2, 3, 4, 6 or 7: input wig file", false, null);
		p.addStringArg("-ot", "For task 1: output file for table of total base counts", false, null);
		p.addStringArg("-on", "For task 2: prefix for separate output wig tracks for each nucleotide", false, null);
		p.addStringArg("-os", "For task 3: output wig file where values have been shifted to a new position in direction of transcription, specified by an offset", false, null);
		p.addIntArg("-of", "For task 3: offset for shifted values (can be negative)", false, -1);
		p.addStringArg("-om", "For task 4: output merged wig file", false, null);
		p.addStringArg("-ob", "For task 5: prefix for output binary wig tracks of transcribed nucleotides", false, null);
		p.addDoubleArg("-mw", "For task 6: minimum wig value to include in intersection", false, 0);
		p.addStringArg("-id", "For task 6: identifier to append to names written to bed file", false, null);
		p.addStringArg("-oi", "For task 6: output bed file of genes intersected with wig data", false, null);
		p.addStringArg("-otd", "For task 7: output file for table of total dinucleotide counts", false, null);
		p.addIntArg("-wfd", "For task 7: first position of dinucleotide relative to wig position", false, 0);
		
		p.parse(args);
		if(p.getFlagsAndValues().isEmpty()) {
			p.printHelpMessage();
			System.exit(-1);
		}
		String genomeFasta = p.getStringArg("-g");
		String bedFile = p.getStringArg("-b");
		String wigFile = p.getStringArg("-w");
		String outTable = p.getStringArg("-ot");
		String outNucTracks = p.getStringArg("-on");
		String outShifted = p.getStringArg("-os");
		int offset = p.getIntArg("-of");
		String wig2 = p.getStringArg("-w2");
		String outMerged = p.getStringArg("-om");
		String outBinary = p.getStringArg("-ob");
		double minWig = p.getDoubleArg("-mw");
		String identifier = p.getStringArg("-id");
		String outIntersected = p.getStringArg("-oi");
		String outDinucTable = p.getStringArg("-otd");
		int dinucFirstPos = p.getIntArg("-wfd");
		
		if(outTable != null || outNucTracks != null || outShifted != null || outDinucTable != null) {
			// Genome and annotation are needed
			if(genomeFasta == null) {
				throw new IllegalArgumentException("Must provide genome with -g option");
			}
			if(bedFile == null) {
				throw new IllegalArgumentException("Must provide annotation with -b option");
			}
			if(wigFile == null) {
				throw new IllegalArgumentException("Must provide wig file with -w option");
			}
			WigTools wsc = new WigTools(wigFile, genomeFasta, bedFile);
			if(outTable != null) wsc.writeBaseCountsToTable(outTable);
			if(outNucTracks != null) wsc.writeTracksByNuceotide(outNucTracks);
			if(outShifted != null) wsc.writeShiftedWig(outShifted, offset);
			if(outDinucTable != null) wsc.writeDinucCountsToTable(dinucFirstPos, outDinucTable);
		}
		
		if(outMerged != null) {
			if(wig2 == null) {
				throw new IllegalArgumentException("Must provide second wig file with -w2 option");
			}
			mergeWigFiles(wigFile, wig2, outMerged);
		}
		
		if(outBinary != null) {
			if(genomeFasta == null) {
				throw new IllegalArgumentException("Must provide genome with -g option");
			}
			if(bedFile == null) {
				throw new IllegalArgumentException("Must provide annotation with -b option");
			}
			WigTools w = new WigTools(genomeFasta, bedFile);
			w.writeBinaryNucleotideTracks(outBinary);
		}
		
		if(outIntersected != null) {
			if(identifier == null) {
				throw new IllegalArgumentException("Must provide identifier with -id option");
			}
			if(bedFile == null) {
				throw new IllegalArgumentException("Must provide annotation with -b option");
			}
			if(wigFile == null) {
				throw new IllegalArgumentException("Must provide wig file with -w option");
			}
			if(genomeFasta == null) {
				throw new IllegalArgumentException("Must provide genome with -g option");
			}
			WigTools w = new WigTools(wigFile, genomeFasta, bedFile);
			w.writeWigOverlappingPositionsAsBed(outIntersected, minWig, identifier);
		}
		
		logger.info("All done.");
		
	}
	
	
	

}
