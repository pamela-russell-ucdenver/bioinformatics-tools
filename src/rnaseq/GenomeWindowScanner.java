/**
 * 
 */
package rnaseq;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.log4j.Logger;

import broad.core.parser.CommandLineParser;
import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import broad.pda.annotation.BEDFileParser;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Annotation.Strand;
import nextgen.core.annotation.Gene;
import nextgen.core.coordinatesystem.GenomicSpace;
import nextgen.core.model.ScanStatisticDataAlignmentModel;
import nextgen.core.model.score.CountScore;
import nextgen.core.model.score.ScanStatisticScore;
import nextgen.core.utils.AlignmentUtils;

/**
 * @author prussell
 *
 */
public class GenomeWindowScanner {

	private ScanStatisticDataAlignmentModel data;
	private GenomicSpace genomeSpace;
	private static double DEFAULT_PVAL_CUTOFF = 0.01;
	private TreeSet<CountScore> windowsOfInterest;
	private static Logger logger = Logger.getLogger(GenomeWindowScanner.class.getName());
	private String bamFile;
	private int windowSize;
	
	private GenomeWindowScanner(String bam, String chrSizeFile, int window, int overlap) throws IOException {
		this(bam, chrSizeFile, window, overlap, DEFAULT_PVAL_CUTOFF, null);
	}
	
	private GenomeWindowScanner(String bam, String chrSizeFile, int window, int overlap, double scanPvalCutoff) throws IOException {
		this(bam, chrSizeFile, window, overlap, scanPvalCutoff, null);
	}	
	
	private GenomeWindowScanner(String bam, String chrSizeFile, int window, int overlap, double scanPvalCutoff, String bedRegionsToScan) throws IOException {
		bamFile = bam;
		windowSize = window;
		genomeSpace = new GenomicSpace(chrSizeFile);
		logger.info("Constructing data alignment model...");
		data = new ScanStatisticDataAlignmentModel(bamFile, genomeSpace, false);
		logger.info("Done constructing data alignment model.");
		windowsOfInterest = new TreeSet<CountScore>();
		if(bedRegionsToScan == null) {
			scanWholeGenome(window, overlap, scanPvalCutoff);
		} else {
			Collection<Gene> regionsToScan = BEDFileParser.loadData(new File(bedRegionsToScan));
			scanRegions(regionsToScan, window, overlap, scanPvalCutoff);
		}
	}
	
	/**
	 * Scan windows in defined regions and keep windows that pass the scan P value cutoff
	 * @param regionsToScan Regions to scan
	 * @param window Window size
	 * @param overlap Overlap
	 * @param scanPvalCutoff Scan P value cutoff
	 */
	private void scanRegions(Collection<? extends Annotation> regionsToScan, int window, int overlap, double scanPvalCutoff) {
		logger.info("Scanning windows across " + regionsToScan.size() + " regions.");
		logger.info("Window size=" + window + ", overlap=" + overlap + ", P-value cutoff=" + scanPvalCutoff);
		int numDone = 0;
		for(Annotation region : regionsToScan) {
			Iterator<CountScore> scoreIter = data.scan(region, window, overlap);
			try {
				while(scoreIter.hasNext()) {
					CountScore countScore = scoreIter.next();
					ScanStatisticScore scanScore = new ScanStatisticScore(data, countScore.getAnnotation(), false);
					numDone++;
					if(numDone % 100000 == 0) {
						logger.info("Scanned " + numDone + " windows of which " + windowsOfInterest.size() + " are significant.");
					}
					double pval = scanScore.getScanPvalue();
					if(pval < scanPvalCutoff) {
						countScore.getAnnotation().setScore(countScore.getCount());
						windowsOfInterest.add(countScore);
					}
					scanScore = null;
				}
			} catch (Exception e) {
				logger.warn("Caught exception, skipping region " + region.getChr() + ":" + region.getStart() + "-" + region.getEnd());
			}
		}
		logger.info("Done scanning windows. There are " + windowsOfInterest.size() + " significant windows.");
	}
	
	/**
	 * Scan all windows in the genome and keep windows that pass the scan P value cutoff
	 * @param window Window size
	 * @param overlap Overlap
	 * @param scanPvalCutoff Scan P value cutoff
	 */
	private void scanWholeGenome(int window, int overlap, double scanPvalCutoff) {
		logger.info("Scanning windows. Window size=" + window + ", overlap=" + overlap + ", P-value cutoff=" + scanPvalCutoff);
		Iterator<ScanStatisticScore> scoreIter = data.scanAll(window, overlap);
		int numDone = 0;
		while(scoreIter.hasNext()) {
			ScanStatisticScore score = scoreIter.next();
			numDone++;
			if(numDone % 1000 == 0) {
				logger.info("Scanned " + numDone + " windows of which " + windowsOfInterest.size() + " are significant.");
			}
			double pval = score.getScanPvalue();
			if(pval < scanPvalCutoff) {
				windowsOfInterest.add(score);
			}
		}
		logger.info("Done scanning windows. There are " + windowsOfInterest.size() + " significant windows.");
	}
	
	/**
	 * Filter windows of interest that overlap an annotation from the set, or windows that don't overlap an annotation from the set
	 * Modifies the set of windows of interest
	 * @param bedFile Bed file of genes to check for overlaps
	 * @param keepOverlappers Whether to keep overlappers (set to true) or keep non-overlappers (set to false)
	 * @throws IOException
	 */
	private void filterWindowsOverlappingAnnotation(String bedFile, boolean keepOverlappers) throws IOException {
		if(!keepOverlappers) logger.info("Removing windows that overlap genes in file " + bedFile + ".");
		if(keepOverlappers) logger.info("Removing windows that do not overlap genes in file " + bedFile + ".");
		Map<String, Collection<Gene>> otherGenes = BEDFileParser.loadDataByChr(new File(bedFile));
		logger.info("Before filtering there are " + windowsOfInterest.size() + " windows.");
		int numDone = 0;
		TreeSet<CountScore> windowsToRemove = new TreeSet<CountScore>();
		for(CountScore score : windowsOfInterest) {
			numDone++;
			if(numDone % 1000 == 0) {
				logger.info("Removed " + windowsToRemove.size() + " of " + numDone + " windows.");
			}
			Gene window = new Gene(score.getAnnotation());
			String chr = window.getReferenceName();
			if(!otherGenes.containsKey(chr)) continue;
			for(Gene otherGene : otherGenes.get(chr)) {
				boolean overlaps = window.overlaps(otherGene, true);
				boolean keep = (overlaps && keepOverlappers) || (!overlaps && !keepOverlappers);
				if(!keep) {
					windowsToRemove.add(score);
				}
			}
		}
		for(CountScore score : windowsToRemove) {
			windowsOfInterest.remove(score);
		}
		logger.info("Done filtering windows. There are now " + windowsOfInterest.size() + " remaining windows.");
	}
	
	/**
	 * Get the windows with most significant scan P values
	 * @param numWindowsToGet Number of top windows to get
	 * @return The most significant windows by scan P value
	 */
	private TreeSet<Annotation> getTopWindows(int numWindowsToGet) {
		logger.info("Getting the top " + numWindowsToGet + " windows...");
		TreeSet<Annotation> rtrn = new TreeSet<Annotation>();
		int numAdded = 0;
		Iterator<CountScore> iter = windowsOfInterest.descendingIterator();
		while(iter.hasNext() && numAdded < numWindowsToGet) {
			CountScore score = iter.next();
			Annotation window = score.getAnnotation();
			rtrn.add(window);
			numAdded++;
		}
		logger.info("Got " + rtrn.size() + " windows.");
		return rtrn;
	}
	
	/**
	 * Write the windows with most significant scan P values to file
	 * @param numWindowsToGet Number of top windows to write
	 * @param outBedFile Output bed file
	 * @throws IOException
	 */
	private void writeTopWindowsBed(int numWindowsToGet, String outBedFile) throws IOException {
		logger.info("Writing the top " + numWindowsToGet + " to bed file " + outBedFile);
		TreeSet<Annotation> windowsToWrite = getTopWindows(numWindowsToGet);
		FileWriter w = new FileWriter(outBedFile);
		for(Annotation window : windowsToWrite) {
			w.write(window.toBED() + "\n");
		}
		w.close();
		logger.info("Done writing top windows to bed file.");
	}
	
	/**
	 * Write fasta sequence file for windows that are already stored in a bed file
	 * @param bedFile The bed file of windows
	 * @param genomeFasta Genome fasta file
	 * @param outFastaFile Output fasta file
	 * @param antisenseSequence Write antisense of each sequence
	 * @param toUpper Convert to upper case
	 * @throws IOException
	 */
	private void writeTopWindowsFastaFromBedFile(String bedFile, String genomeFasta, String outFastaFile, boolean antisenseSequence, boolean toUpper) throws IOException {
		Collection<Gene> genes = BEDFileParser.loadData(new File(bedFile));
		logger.info("Writing " + genes.size() + " windows from file " + bedFile + " to fasta file " + outFastaFile);
		TreeSet<Annotation> windows = new TreeSet<Annotation>();
		for(Gene gene : genes) {
			windows.add(gene);
		}
		writeTopWindowsFasta(windows, genomeFasta, outFastaFile, antisenseSequence, toUpper);
	}
	
	/**
	 * Write fasta sequence file for top scoring windows
	 * @param genomeFasta Genome fasta file
	 * @param numWindowsToGet Number of top windows to write
	 * @param outFastaFile Output fasta file
	 * @param antisenseSequence Write antisense of each sequence
	 * @param toUpper Convert to upper case
	 * @throws IOException
	 */
	private void writeTopWindowsFasta(String genomeFasta, int numWindowsToGet, String outFastaFile, boolean antisenseSequence, boolean toUpper) throws IOException {
		logger.info("Writing the top " + numWindowsToGet + " to fasta file " + outFastaFile);
		TreeSet<Annotation> windowsToWrite = getTopWindows(numWindowsToGet);
		writeTopWindowsFasta(windowsToWrite, genomeFasta, outFastaFile, antisenseSequence, toUpper);
	}
	
	/**
	 * Write fasta sequence file for set of windows
	 * @param windowsToWrite The windows to write
	 * @param genomeFasta Genome fasta file
	 * @param outFastaFile Output fasta file
	 * @param antisenseSequence Write antisense of each sequence
	 * @param toUpper Convert to upper case
	 * @throws IOException
	 */
	private void writeTopWindowsFasta(TreeSet<Annotation> windowsToWrite, String genomeFasta, String outFastaFile, boolean antisenseSequence, boolean toUpper) throws IOException {
		
		// Load chromosomes
		logger.info("Loading chromosomes from file " + genomeFasta + "...");
		FastaSequenceIO fsio = new FastaSequenceIO(genomeFasta);
		List<Sequence> chrList = fsio.loadAll();
		Map<String, Sequence> chromosomes = new TreeMap<String, Sequence>();
		for(Sequence chr : chrList) {
			chromosomes.put(chr.getId(), chr);
		}
		logger.info("Loaded " + chromosomes.size() + " chromosomes.");
		
		// Construct sequences
		logger.info("Constructing " + windowsToWrite.size() + " fasta sequences...");
		int numDone = 0;
		List<Sequence> seqsToWrite = new ArrayList<Sequence>();
		for(Annotation window : windowsToWrite) {
			numDone++;
			if(numDone % 1000 == 0) {
				logger.info("Constructed " + numDone + " sequences.");
			}
			Gene gene = new Gene(window);
			Sequence seq = chromosomes.get(gene.getChr()).getSubsequence(gene);
			boolean minusStrand = window.getOrientation().equals(Strand.NEGATIVE);
			if(antisenseSequence) seq.reverse();
			if(toUpper) seq.uppercase();
			String name = gene.getChr() + ":" + gene.getStart() + "-" + gene.getEnd() + ":";
			if(minusStrand) {
				name += "-";
			} else {
				name += "+";
			}
			seq.setId(name);
			seqsToWrite.add(seq);
		}
		
		// Write fasta file
		FileWriter w = new FileWriter(outFastaFile);
		BufferedWriter b = new BufferedWriter(w);
		FastaSequenceIO fastaWriter = new FastaSequenceIO();
		fastaWriter.write(seqsToWrite, b, windowSize);
		b.close();
		logger.info("Done writing top windows to fasta file.");
	}
	
	/**
	 * Write all significant windows to file
	 * @param outBedFile Output bed file
	 * @throws IOException
	 */
	private void writeAllSignificantWindows(String outBedFile) throws IOException {
		logger.info("Writing significant windows to file " + outBedFile);
		FileWriter w = new FileWriter(outBedFile);
		for(CountScore score : windowsOfInterest) {
			Annotation window = score.getAnnotation();
			w.write(window.toBED() + "\n");
		}
		w.close();
		logger.info("Done writing  windows to file.");		
	}
	
	/**
	 * Assign orientation to each window based on number of reads mapping in each orientation
	 */
	private void assignOrientationToWindows(boolean firstReadTranscriptionStrand, double cutoff) {
		logger.info("Assigning orientation to " + windowsOfInterest.size() + " windows.");
		int numDone = 0;
		for(CountScore score : windowsOfInterest) {
			Annotation window = score.getAnnotation();
			window.setOrientation(AlignmentUtils.assignOrientationToWindow(bamFile, window, firstReadTranscriptionStrand, cutoff));
			if(numDone % 1000 == 0) {
				logger.info("Assigned orientation to " + numDone + " windows.");
			}
		}
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {

		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-b", "Bam file", true);
		p.addStringArg("-s", "Chromosome size file", true);
		p.addIntArg("-w", "Window size", true);
		p.addIntArg("-o", "Window overlap", true);
		p.addDoubleArg("-p", "Scan P value cutoff for window", false, DEFAULT_PVAL_CUTOFF);
		p.addStringArg("-gk", "Bed file of genes to filter against: only keep windows that overlap a gene. Cannot be used in conjunction with -gd.", false);
		p.addStringArg("-gd", "Bed file of genes to filter against: only keep windows that do not overlap a gene. Cannot be used in conjunction with -gk.", false);
		p.addStringArg("-to", "Write top windows to this bed file. Must provide -t.", false);
		p.addIntArg("-t", "The number of top windows to write to bed file", false);
		p.addStringArg("-ao", "Write all significant windows to this bed file", false);
		p.addStringArg("-r", "Bed file of regions to scan", false, null);
		p.addStringArg("-fs", "Write top windows as fasta sequences to this file. Sense direction to orientation of window. Requires -g and -t. Cannot be used in conjuction with -fa", false);
		p.addStringArg("-fa", "Write top windows as fasta sequences to this file. Antisense direction to orientation of window. Requires -g and -t. Cannot be used in conjuction with -fs", false);
		p.addStringArg("-bw", "Bed file of windows to write with -ba or -bs", false);
		p.addStringArg("-ba", "Write windows stored in bed file as fasta sequences to this file. Antisense direction to orientation of window. Requires -g and -bw. Cannot be used in conjunction with -bs.", false);
		p.addStringArg("-bs", "Write windows stored in bed file as fasta sequences to this file. Sense direction to orientation of window. Requires -g and -bw. Cannot be used in conjunction with -ba.", false);
		p.addBooleanArg("-st", "First read is transcription strand", true);
		p.addStringArg("-g", "Genome fasta. Required for -fs and -fa.", false);
		p.addBooleanArg("-u", "Write fasta sequences in all uppercase", false, true);
		p.addDoubleArg("-c", "Cutoff for proportion of reads to assign orientation to window", false, 0.5);
		p.parse(args);
		
		// Check for invalid arguments
		if(p.hasStringFlag("-gk") && p.hasStringFlag("-gd")) {
			throw new IllegalArgumentException("Cannot specify both -gk and -gd.");
		}
		if(p.hasStringFlag("-fs") && p.hasStringFlag("-fa")) {
			throw new IllegalArgumentException("Cannot specify both -fs and -fa.");
		}
		if(p.hasStringFlag("-bs") && p.hasStringFlag("-ba")) {
			throw new IllegalArgumentException("Cannot specify both -bs and -ba.");
		}
		if(p.hasStringFlag("-to") && !p.hasIntFlag("-t")) {
			throw new IllegalArgumentException("If using -to must also specify -t.");
		}
		if(p.hasStringFlag("-fs") && !p.hasIntFlag("-t")) {
			throw new IllegalArgumentException("If using -fs must also specify -t.");
		}
		if(p.hasStringFlag("-fa") && !p.hasIntFlag("-t")) {
			throw new IllegalArgumentException("If using -fa must also specify -t.");
		}
		if(p.hasStringFlag("-fs") && !p.hasStringFlag("-g")) {
			throw new IllegalArgumentException("If using -fs must also specify -g.");
		}
		if(p.hasStringFlag("-fa") && !p.hasStringFlag("-g")) {
			throw new IllegalArgumentException("If using -fa must also specify -g.");
		}
		if(p.hasStringFlag("-bs") && !p.hasStringFlag("-g")) {
			throw new IllegalArgumentException("If using -bs must also specify -g.");
		}
		if(p.hasStringFlag("-ba") && !p.hasStringFlag("-g")) {
			throw new IllegalArgumentException("If using -ba must also specify -g.");
		}
		if(p.hasStringFlag("-bs") && !p.hasStringFlag("-bw")) {
			throw new IllegalArgumentException("If using -bs must also specify -bw.");
		}
		if(p.hasStringFlag("-ba") && !p.hasStringFlag("-bw")) {
			throw new IllegalArgumentException("If using -ba must also specify -bw.");
		}
		
		// Construct object
		String bamFile = p.getStringArg("-b");
		String chrSizeFile = p.getStringArg("-s");
		int windowSize = p.getIntArg("-w");
		int overlap = p.getIntArg("-o");
		double scanPvalCutoff = p.getDoubleArg("-p");
		String regionsBedFile = p.getStringArg("-r");
		double cutoff = p.getDoubleArg("-c");
		GenomeWindowScanner gws = new GenomeWindowScanner(bamFile, chrSizeFile, windowSize, overlap, scanPvalCutoff, regionsBedFile);
		
		// Filter for gene overlap
		if(p.hasStringFlag("-gk")) {
			String bedFile = p.getStringArg("-gk");
			gws.filterWindowsOverlappingAnnotation(bedFile, true);
		}
		if(p.hasStringFlag("-gd")) {
			String bedFile = p.getStringArg("-gd");
			gws.filterWindowsOverlappingAnnotation(bedFile, false);
		}
		
		// Assign orientation to each window
		boolean firstReadTranscriptionStrand = p.getBooleanArg("-st");
		gws.assignOrientationToWindows(firstReadTranscriptionStrand, cutoff);
		
		// Write all significant windows to file
		if(p.hasStringFlag("-ao")) {
			String bedFile = p.getStringArg("-ao");
			gws.writeAllSignificantWindows(bedFile);
		}
		
		// Write top significant windows to file
		if(p.hasStringFlag("-to")) {
			String bedFile = p.getStringArg("-to");
			int numToGet = p.getIntArg("-t");
			gws.writeTopWindowsBed(numToGet, bedFile);
		}
		
		// Write top significant windows as fasta
		if(p.hasStringFlag("-fs")) {
			boolean toUpper = p.getBooleanArg("-u");
			String fastaFile = p.getStringArg("-fs");
			String genomeFasta = p.getStringArg("-g");
			int numToGet = p.getIntArg("-t");
			gws.writeTopWindowsFasta(genomeFasta, numToGet, fastaFile, false, toUpper);
		}
		if(p.hasStringFlag("-fa")) {
			boolean toUpper = p.getBooleanArg("-u");
			String fastaFile = p.getStringArg("-fa");
			String genomeFasta = p.getStringArg("-g");
			int numToGet = p.getIntArg("-t");
			gws.writeTopWindowsFasta(genomeFasta, numToGet, fastaFile, true, toUpper);
		}
		
		// Write windows from bed file as fasta
		if(p.hasStringFlag("-bs")) {
			boolean toUpper = p.getBooleanArg("-u");
			String bedFile = p.getStringArg("-bw");
			String fastaFile = p.getStringArg("-bs");
			String genomeFasta = p.getStringArg("-g");
			gws.writeTopWindowsFastaFromBedFile(bedFile, genomeFasta, fastaFile, false, toUpper);
		}
		if(p.hasStringFlag("-ba")) {
			boolean toUpper = p.getBooleanArg("-u");
			String bedFile = p.getStringArg("-bw");
			String fastaFile = p.getStringArg("-ba");
			String genomeFasta = p.getStringArg("-g");
			gws.writeTopWindowsFastaFromBedFile(bedFile, genomeFasta, fastaFile, true, toUpper);
		}
		
	}

}
