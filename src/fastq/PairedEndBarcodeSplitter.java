package fastq;

import guttmanlab.core.alignment.SmithWatermanAlignment;
import guttmanlab.core.util.CommandLineParser;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;

import nextgen.core.utils.FileUtil;

import org.apache.log4j.Logger;

import broad.pda.seq.fastq.FastqParser;
import broad.pda.seq.fastq.FastqSequence;

/**
 * Barcode splitter for paired end reads
 * Keeps the fastq files in same order
 * @author prussell
 *
 */
public class PairedEndBarcodeSplitter {
	
	private String read1fq; // Original read1 fastq file
	private String read2fq; // Original read2 fastq file
	private Map<String, String> read1out; // Output files for read1 by barcode
	private Map<String, String> read2out; // Output files for read2 by barcode
	private Map<String, BufferedWriter> read1writers; // Read1 writers by barcode
	private Map<String, BufferedWriter> read2writers; // Read2 writers by barcode
	private BufferedWriter unmatched1writer; // Writer for unmatched read1s
	private BufferedWriter unmatched2writer; // Writer for unmatched read2s
	private Collection<String> barcodes; // All possible barcodes
	private boolean removeBarcode; // Whether to remove the barcode from reads before writing
	private FastqParser fastq1parser; // Read1 fastq reader
	private FastqParser fastq2parser; // Read2 fastq reader
	private BarcodeMatcher matcher; // Barcode matcher implementation
	
	private static Logger logger = Logger.getLogger(PairedEndBarcodeSplitter.class.getName());

	/**
	 * Mate 1 or 2, both or neither
	 * @author prussell
	 *
	 */
	private enum Mate {
		MATE1,
		MATE2,
		BOTH,
		NEITHER;
		
		public String toString() {
			switch(this) {
			case BOTH:
				return "BOTH";
			case MATE1:
				return "MATE1";
			case MATE2:
				return "MATE2";
			case NEITHER:
				return "NEITHER";
			default:
				throw new IllegalArgumentException("Not supported");
			}
		}
		
	}
	
	/**
	 * @param barcodeMatcher Barcode matcher implementation
	 * @param read1fastq Original read 1 fastq file
	 * @param read2fastq Original read 2 fastq file
	 * @param removeBarcodeFromReads Remove barcode from reads before writing
	 * @throws IOException
	 */
	private PairedEndBarcodeSplitter(BarcodeMatcher barcodeMatcher, String read1fastq, String read2fastq, boolean removeBarcodeFromReads) throws IOException {
		matcher = barcodeMatcher;
		matcher.initializeBarcodes(this);
		read1fq = read1fastq;
		read2fq = read2fastq;
		removeBarcode = removeBarcodeFromReads;
		resetReaders();
		resetWriters();
	}
	
	public PairedEndBarcodeSplitter() {}

	/**
	 * Reset all the file writers
	 * @throws IOException
	 */
	private void resetWriters() throws IOException {
		logger.info("Resetting writers");
		read1out = new HashMap<String, String>();
		read2out = new HashMap<String, String>();
		read1writers = new HashMap<String, BufferedWriter>();
		read2writers = new HashMap<String, BufferedWriter>();
		for(String barcode : barcodes) {
			read1out.put(barcode, addBarcodeToFileName(read1fq, barcode));
			read2out.put(barcode, addBarcodeToFileName(read2fq, barcode));
			read1writers.put(barcode, new BufferedWriter(new FileWriter(read1out.get(barcode))));
			read2writers.put(barcode, new BufferedWriter(new FileWriter(read2out.get(barcode))));
			logger.info("Writing reads with barcode " + barcode + " to " + read1out.get(barcode) + " and " + read2out.get(barcode));
		}
		unmatched1writer = new BufferedWriter(new FileWriter("unmatched_1.fq"));
		unmatched2writer = new BufferedWriter(new FileWriter("unmatched_2.fq"));
	}
	
	/**
	 * Close all the file writers
	 * @throws IOException
	 */
	private void closeWriters() throws IOException {
		logger.info("Closing writers");
		for(BufferedWriter writer : read1writers.values()) {
			writer.close();
		}
		for(BufferedWriter writer : read2writers.values()) {
			writer.close();
		}
		unmatched1writer.close();
		unmatched2writer.close();
	}
	
	/**
	 * Reset the fastq parsers
	 * @throws IOException 
	 */
	private void resetReaders() throws IOException {
		logger.info("Resetting readers");
		fastq1parser = new FastqParser();
		fastq1parser.start(new File(read1fq));
		fastq2parser = new FastqParser();
		fastq2parser.start(new File(read2fq));
	}
	
	/**
	 * Add barcode to output file name, moving the file extension to the end
	 * @param fileName File name to add barcode to
	 * @param barcode Barcode to add
	 * @return File name with barcode added and extension moved
	 */
	private static String addBarcodeToFileName(String fileName, String barcode) {
		if(fileName.contains(".fq")) {
			return fileName.replaceAll(".fq", "." + barcode + ".fq");
		}
		if(fileName.contains(".fastq")) {
			return fileName.replaceAll(".fastq", "." + barcode + ".fastq");
		}
		return fileName + "." + barcode;
	}
	
	/**
	 * Write matches to separate fastq files by barcode
	 * @throws IOException
	 */
	private void writeMatches() throws IOException {
		logger.info("");
		logger.info("Writing separate files by barcode...");
		resetReaders();
		resetWriters();
		int numDone = 0;
		while(fastq1parser.hasNext() && fastq2parser.hasNext()) {
			numDone++;
			if(numDone % 1000000 == 0) {
				logger.info("Finished " + numDone + " records.");
			}
			FastqSequence read1 = fastq1parser.next();
			FastqSequence read2 = fastq2parser.next();
			String barcode = matcher.identifyBarcode(read1, read2);
			if(barcode != null) {
				FastqSequence toWrite1 = matcher.getRead1SeqToWrite(read1, barcode);
				FastqSequence toWrite2 = matcher.getRead2SeqToWrite(read2, barcode);
				toWrite1.write(read1writers.get(barcode));
				toWrite2.write(read2writers.get(barcode));
			} else {
				read1.write(unmatched1writer);
				read2.write(unmatched2writer);
			}
		}
		closeWriters();
		logger.info("Done writing files.");
	}
	
	/**
	 * Remove barcode from fastq record
	 * @param record Record
	 * @param barcode Barcode to remove
	 * @return A new fastq record with the barcode removed
	 */
	private static FastqSequence removeBarcode(FastqSequence record, String barcode) {
		throw new UnsupportedOperationException("not implemented");
	}
	
	/**
	 * A way to match barcodes to paired reads
	 * @author prussell
	 *
	 */
	private interface BarcodeMatcher {
		
		/**
		 * Whether the read contains the barcode
		 * @param read Read
		 * @param barcode Barcode
		 * @return True iff the read contains the barcode
		 */
		public boolean readHasBarcode(FastqSequence read, String barcode);
		
		/**
		 * Get the barcode for the read pair
		 * @param read1 Read 1
		 * @param read2 Read 2
		 * @return The barcode for the read pair or null if no match
		 */

		public String identifyBarcode(FastqSequence read1, FastqSequence read2);
		
		/**
		 * Initialize the set of barcodes for the PairedEndBarcodeSplitter object
		 * @throws IOException 
		 */
		public void initializeBarcodes(PairedEndBarcodeSplitter barcodeSplitter) throws IOException;
		
		/**
		 * Get read1 sequence to write to output
		 * @param read1 Read1
		 * @param barcode A barcode contained in read1
		 * @return Read1 sequence to write
		 */
		public FastqSequence getRead1SeqToWrite(FastqSequence read1, String barcode);
		
		/**
		 * Get read2 sequence to write to output
		 * @param read2 Read2
		 * @param barcode A barcode contained in read1
		 * @return Read2 sequence to write
		 */
		public FastqSequence getRead2SeqToWrite(FastqSequence read2, String barcode);
		
	}
	
	/**
	 * Barcode matcher where either read1 or read2 contains a mate-specific barcode, but not both
	 * E.g. RNA and DNA have separate barcode ligations on different ends of the molecule;
	 * DPM sequence appears in read1 and RPM sequence appears in read2; no read pair should have both
	 * @author prussell
	 *
	 */
	private class OneMateContainsBarcode implements BarcodeMatcher {
		
		private Map<Mate, String> mateSpecificBarcode;
		private int maxMismatchRead1barcode;
		private int maxMismatchRead2barcode;
		
		/**
		 * @param mate1barcode Barcode that can be contained in read1
		 * @param mate2barcode Barcode that can be contained in read2
		 * @param maxMismatchMate1barcode Max mismatch in read1 barcode
		 * @param maxMismatchMate2barcode Max mismatch in read2 barcode
		 */
		public OneMateContainsBarcode(String mate1barcode, String mate2barcode, int maxMismatchMate1barcode, int maxMismatchMate2barcode) {
			mateSpecificBarcode = new HashMap<Mate, String>();
			mateSpecificBarcode.put(Mate.MATE1, mate1barcode);
			mateSpecificBarcode.put(Mate.MATE2, mate2barcode);
			maxMismatchRead1barcode = maxMismatchMate1barcode;
			maxMismatchRead2barcode = maxMismatchMate2barcode;
		}
		
		@Override
		public boolean readHasBarcode(FastqSequence read, String barcode) {
			throw new UnsupportedOperationException();
		}

		/**
		 * Whether the read contains the barcode
		 * @param read Read
		 * @param barcode Barcode
		 * @param maxMismatches Max mismatches
		 * @return True iff the read contains the barcode with at most the max mismatches
		 */
		private boolean readHasBarcode(FastqSequence read, String barcode, int maxMismatches) {
			return SmithWatermanAlignment.containsFullLengthUngappedMatch(read.getSequence(), barcode, maxMismatches);
		}
		
		@Override
		public String identifyBarcode(FastqSequence read1, FastqSequence read2) {
			boolean read1hasSpecificBarcode = readHasBarcode(read1, mateSpecificBarcode.get(Mate.MATE1), maxMismatchRead1barcode);
			boolean read2hasSpecificBarcode = readHasBarcode(read2, mateSpecificBarcode.get(Mate.MATE2), maxMismatchRead2barcode);
			if(read1hasSpecificBarcode) {
				if(read2hasSpecificBarcode) {
					return Mate.BOTH.toString();
				}
				return mateSpecificBarcode.get(Mate.MATE1);
			}
			if(read2hasSpecificBarcode) {
				return mateSpecificBarcode.get(Mate.MATE2);
			}
			return Mate.NEITHER.toString();
		}

		@Override
		public void initializeBarcodes(PairedEndBarcodeSplitter barcodeSplitter) throws IOException {
			barcodeSplitter.barcodes = new HashSet<String>();
			barcodeSplitter.barcodes.addAll(mateSpecificBarcode.values());
			barcodeSplitter.barcodes.add(Mate.BOTH.toString());
			barcodeSplitter.barcodes.add(Mate.NEITHER.toString());
		}

		@Override
		public FastqSequence getRead1SeqToWrite(FastqSequence read1, String barcode) {
			return read1;
		}

		@Override
		public FastqSequence getRead2SeqToWrite(FastqSequence read2, String barcode) {
			return read2;
		}
		
	}
		
	/**
	 * Check whether a fastq record contains a perfect match of a sequence
	 * @author prussell
	 *
	 */
	private class PerfectMatch implements BarcodeMatcher {
		
		private String barcodeListFile;
		private boolean read1hasBarcode;

		/**
		 * @param barcodeFile File containing simple list of barcodes
		 * @param read1Barcode True iff read1 has the barcode
		 */
		public PerfectMatch(String barcodeFile, boolean read1Barcode) {
			barcodeListFile = barcodeFile;
			read1hasBarcode = read1Barcode;
		}
		
		@Override
		public boolean readHasBarcode(FastqSequence read, String barcode) {
			String readSeq = read.getSequence();
			return readSeq.contains(barcode);
		}

		@Override
		public void initializeBarcodes(PairedEndBarcodeSplitter barcodeSplitter) throws IOException {
			barcodeSplitter.barcodes = new HashSet<String>();
			barcodeSplitter.barcodes.addAll(FileUtil.fileLinesAsList(barcodeListFile));
		}
	
		@Override
		public String identifyBarcode(FastqSequence read1, FastqSequence read2) {
			FastqSequence barcodeRead = read1hasBarcode ? read1 : read2;
			for(String barcode : barcodes) {
				if(readHasBarcode(barcodeRead, barcode)) {
					return barcode;
				}
			}
			return null;
		}

		@Override
		public FastqSequence getRead1SeqToWrite(FastqSequence read1, String barcode) {
			return (read1hasBarcode && removeBarcode) ? removeBarcode(read1, barcode) : read1;
		}

		@Override
		public FastqSequence getRead2SeqToWrite(FastqSequence read2, String barcode) {
			return (!read1hasBarcode && removeBarcode) ? removeBarcode(read2, barcode) : read2;
		}

	}
	
	/**
	 * The options for barcode matcher implementations
	 * @author prussell
	 *
	 */
	private enum MatchImplementation {
		
		PERFECT_MATCH,
		ONE_MATE_CONTAINS_BARCODE;
		
		public String toString() {
			switch(this) {
			case ONE_MATE_CONTAINS_BARCODE:
				return "one_mate_contains_barcode";
			case PERFECT_MATCH:
				return "perfect_match";
			default:
				throw new IllegalArgumentException("Not implemented");
			}
		}
		
		public static MatchImplementation fromString(String s) {
			if(s.equals(PERFECT_MATCH.toString())) {
				return PERFECT_MATCH;
			}
			if(s.equals(ONE_MATE_CONTAINS_BARCODE.toString())) {
				return ONE_MATE_CONTAINS_BARCODE;
			}
			throw new IllegalArgumentException("Not implemented");
		}
		
		public static String commaSeparatedList() {
			String rtrn = "";
			for(int i = 0; i < MatchImplementation.values().length; i++) {
				if(i > 0) {
					rtrn += ",";
				}
				rtrn += MatchImplementation.values()[i].toString();
			}
			return rtrn;
		}
		
	}
	
	/**
	 * Get barcode matcher object specified on the command line
	 * @param p Command line parser
	 * @return Barcode matcher object
	 */
	private static BarcodeMatcher getBarcodeMatcher(CommandLineParser p) {
		String barcodeFile = p.getStringArg(perfectMatchBarcodeFileOption);
		String barcodeImpl = p.getStringArg(barcodeMatcherOption);
		boolean perfectMatchRead1hasBarcode = p.getBooleanArg(perfectMatchBarcodeOnRead1Option);
		String mate1barcode = p.getStringArg(oneMateContainsBarcodeMate1BarcodeOption);
		String mate2barcode = p.getStringArg(oneMateContainsBarcodeMate2BarcodeOption);
		int maxMismatchMate1barcode = p.getIntArg(oneMateContainsBarcodeMaxMismatchRead1Option);
		int maxMismatchMate2barcode = p.getIntArg(oneMateContainsBarcodeMaxMismatchRead2Option);
		
		MatchImplementation implementation = MatchImplementation.fromString(barcodeImpl);
		
		switch(implementation) {
		case ONE_MATE_CONTAINS_BARCODE:
			if(mate1barcode == null) {
				throw new IllegalArgumentException("Must provide " + oneMateContainsBarcodeMate1BarcodeOption);
			}
			if(mate2barcode == null) {
				throw new IllegalArgumentException("Must provide " + oneMateContainsBarcodeMate2BarcodeOption);
			}
			if(maxMismatchMate1barcode < 0) {
				throw new IllegalArgumentException("Must provide " + oneMateContainsBarcodeMaxMismatchRead1Option);
			}
			if(maxMismatchMate2barcode < 0) {
				throw new IllegalArgumentException("Must provide " + oneMateContainsBarcodeMaxMismatchRead2Option);
			}
			return new PairedEndBarcodeSplitter().new OneMateContainsBarcode(mate1barcode, mate2barcode, maxMismatchMate1barcode, maxMismatchMate2barcode);
		case PERFECT_MATCH:
			if(barcodeFile == null) {
				throw new IllegalArgumentException("Must provide " + perfectMatchBarcodeFileOption);
			}
			return new PairedEndBarcodeSplitter().new PerfectMatch(barcodeFile, perfectMatchRead1hasBarcode);
		default:
			throw new IllegalArgumentException("Not implemented: " + implementation.toString());
		}
		
	}
	
	private static String perfectMatchBarcodeFileOption = "-pmb";
	private static String read1fqOption = "-f1";
	private static String read2fqOption = "-f2";
	private static String perfectMatchBarcodeOnRead1Option = "-pmb1";
	private static String removeBarcodeOption = "-rb";
	private static String barcodeMatcherOption = "-bmi";
	private static String oneMateContainsBarcodeMate1BarcodeOption = "-omb1";
	private static String oneMateContainsBarcodeMate2BarcodeOption = "-omb2";
	private static String oneMateContainsBarcodeMaxMismatchRead1Option = "-omm1";
	private static String oneMateContainsBarcodeMaxMismatchRead2Option = "-omm2";
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg(perfectMatchBarcodeFileOption, "For perfect match implementation, file containing simple list of barcodes", false, null);
		p.addStringArg(read1fqOption, "Read1 fastq file", true);
		p.addStringArg(read2fqOption, "Read2 fastq file", true);
		p.addBooleanArg(perfectMatchBarcodeOnRead1Option, "For perfect match implementation, true if barcode is on read1, false if read2", false, true);
		p.addBooleanArg(removeBarcodeOption, "Remove barcode when writing separated files", false, true);
		p.addStringArg(barcodeMatcherOption, "Barcode matcher implementation. Options: " + MatchImplementation.commaSeparatedList(), true);
		p.addStringArg(oneMateContainsBarcodeMate1BarcodeOption, "For one mate implementation, read1 barcode", false, null);
		p.addStringArg(oneMateContainsBarcodeMate2BarcodeOption, "For one mate implementation, read2 barcode", false, null);
		p.addIntArg(oneMateContainsBarcodeMaxMismatchRead1Option, "For one mate implementation, max mismatches in read 1 barcode", false, -1);
		p.addIntArg(oneMateContainsBarcodeMaxMismatchRead2Option, "For one mate implementation, max mismatches in read 2 barcode", false, -1);
		p.parse(args);
		String read1fastq = p.getStringArg("-f1");
		String read2fastq = p.getStringArg("-f2");
		boolean removeBarcodeFromReads = p.getBooleanArg("-rb");
		
		BarcodeMatcher matcher = getBarcodeMatcher(p);
		PairedEndBarcodeSplitter bs = new PairedEndBarcodeSplitter(matcher, read1fastq, read2fastq, removeBarcodeFromReads);
		bs.writeMatches();
		
		logger.info("");
		logger.info("All done.");

	}

}
