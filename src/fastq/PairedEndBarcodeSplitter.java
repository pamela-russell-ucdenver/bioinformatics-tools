package fastq;

import guttmanlab.core.alignment.LocalAlignmentUtils;
import guttmanlab.core.alignment.SmithWatermanAlignment;
import guttmanlab.core.util.CommandLineParser;
import guttmanlab.core.util.StringParser;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import nextgen.core.utils.FileUtil;

import org.apache.log4j.Logger;

import broad.core.sequence.Sequence;
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
	private boolean trimReads; // Whether to remove the barcode from reads before writing
	private FastqParser fastq1parser; // Read1 fastq reader
	private FastqParser fastq2parser; // Read2 fastq reader
	private BarcodeMatcher matcher; // Barcode matcher implementation
	
	@SuppressWarnings("unused")
	private static Character[] alphabet = {'A', 'C', 'G', 'T', 'N'};
	
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
		
		/**
		 * Get instance based on string name
		 * @param str Name
		 * @return Instance corresponding to name
		 */
		public static Mate fromString(String str) {
			if(str.equals(MATE1.toString())) {
				return MATE1;
			}
			if(str.equals(MATE2.toString())) {
				return MATE2;
			}
			if(str.equals(BOTH.toString())) {
				return BOTH;
			}
			if(str.equals(NEITHER.toString())) {
				return NEITHER;
			}
			throw new IllegalArgumentException("Not implemented for " + str);
		}
		
	}
	
	/**
	 * @param barcodeMatcher Barcode matcher implementation
	 * @param read1fastq Original read 1 fastq file
	 * @param read2fastq Original read 2 fastq file
	 * @param trimBarcodeFromReads Trim reads before writing
	 * @throws IOException
	 */
	private PairedEndBarcodeSplitter(BarcodeMatcher barcodeMatcher, String read1fastq, String read2fastq, boolean trimBarcodeFromReads) throws IOException {
		matcher = barcodeMatcher;
		matcher.initializeBarcodes(this);
		read1fq = read1fastq;
		read2fq = read2fastq;
		trimReads = trimBarcodeFromReads;
		resetReaders();
		resetWriters();
	}
	
	private PairedEndBarcodeSplitter() {}

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
		
		/**
		 * Trim barcodes, etc. from read1 and get trimmed read
		 * @param read1 Full read1
		 * @param barcode Barcode identified in read1
		 * @return Trimmed version of read1
		 */
		public FastqSequence trimRead1(FastqSequence read1, String barcode);
		
		/**
		 * Trim barcodes, etc. from read2 and get trimmed read
		 * @param read2 Full read2
		 * @param barcode Barcode identified in read2
		 * @return Trimmed version of read2
		 */
		public FastqSequence trimRead2(FastqSequence read2, String barcode);
		
	}
	
	private static Map<String, Integer> posMapFromFile(String file) throws IOException {
		BufferedReader b = new BufferedReader(new FileReader(file));
		StringParser s = new StringParser();
		Map<String, Integer> rtrn = new HashMap<String, Integer>();
		while(b.ready()) {
			s.parse(b.readLine());
			if(s.getFieldCount() != 2) {
				b.close();
				throw new IllegalArgumentException("Line format: string   int");
			}
			rtrn.put(s.asString(0), Integer.valueOf(s.asInt(1)));
		}
		b.close();
		return rtrn;
	}

	
	/**
	 * All barcodes are in one mate or the other
	 * Each barcode has a known position where it appears in the reads
	 * Different barcodes can have different positions
	 * @author prussell
	 *
	 */
	private class OneMateContainsBarcodesAtKnownPositions implements BarcodeMatcher {
		
		private Mate mateWithBarcode;
		private Map<String, Integer> barcodePositions;
		private int maxMismatches;
		
		public OneMateContainsBarcodesAtKnownPositions(Mate mateContainingBarcode, String barcodeFile, int maxMismatch) throws IOException {
			this(mateContainingBarcode, posMapFromFile(barcodeFile), maxMismatch);
		}
		
		public OneMateContainsBarcodesAtKnownPositions(Mate mateContainingBarcode, Map<String, Integer> barcodePos, int maxMismatch) {
			if(mateContainingBarcode != Mate.MATE1 && mateContainingBarcode != Mate.MATE2) {
				throw new IllegalArgumentException("Mate must be mate1 or mate2");
			}
			mateWithBarcode = mateContainingBarcode;
			barcodePositions = barcodePos;
			maxMismatches = maxMismatch;
		}
		
		
		@Override
		public boolean readHasBarcode(FastqSequence read, String barcode) {
			throw new IllegalArgumentException("Not implemented because can't identify whether read is mate 1 or 2");
		}

		@Override
		public String identifyBarcode(FastqSequence read1, FastqSequence read2) {
			String seq = mateWithBarcode == Mate.MATE1 ? read1.getSequence() : read2.getSequence();
			boolean found = false;
			String[] rtrn = new String[1];
			for(String barcode : barcodePositions.keySet()) {
				if(LocalAlignmentUtils.containsFullLengthUngappedMatchAtPosition(seq, barcode, barcodePositions.get(seq).intValue(), maxMismatches, true)) {
					if(found) {
						throw new IllegalStateException("Multiple barcodes match sequence " + seq + " at position " + barcodePositions.get(seq).toString() + " with " + maxMismatches + " allowable mismatches");
					}
					rtrn[0] = barcode;
				}
			}
			return rtrn[0];
		}

		@Override
		public void initializeBarcodes(PairedEndBarcodeSplitter barcodeSplitter) throws IOException {
			throw new UnsupportedOperationException("Not implemented");
		}

		@Override
		public FastqSequence getRead1SeqToWrite(FastqSequence read1, String barcode) {
			if(mateWithBarcode == Mate.MATE2) {
				return read1;
			}
			return trimReads ? trimRead1(read1, barcode) : read1;
		}

		@Override
		public FastqSequence getRead2SeqToWrite(FastqSequence read2, String barcode) {
			if(mateWithBarcode == Mate.MATE1) {
				return read2;
			}
			return trimReads ? trimRead2(read2, barcode) : read2;
		}

		@Override
		public FastqSequence trimRead1(FastqSequence read1, String barcode) {
			if(mateWithBarcode == Mate.MATE2) {
				return read1;
			}
			return getSequenceAfterUngappedMatch(read1, barcode, barcodePositions.get(barcode), maxMismatches);
		}

		@Override
		public FastqSequence trimRead2(FastqSequence read2, String barcode) {
			if(mateWithBarcode == Mate.MATE1) {
				return read2;
			}
			return getSequenceAfterUngappedMatch(read2, barcode, barcodePositions.get(barcode), maxMismatches);
		}
		
	}
	
	/**
	 * Barcode matcher where either read1 or read2 contains a mate-specific barcode, but not both
	 * E.g. RNA and DNA have separate barcode ligations on different ends of the molecule;
	 * DPM sequence appears in read1 and RPM sequence appears in read2; no read pair should have both
	 * @author prussell
	 *
	 */
	@SuppressWarnings("unused")
	private class EitherMateContainsBarcode implements BarcodeMatcher {
		
		protected Map<Mate, String> mateSpecificBarcode;
		protected int maxMismatchRead1barcode;
		protected int maxMismatchRead2barcode;
		protected boolean read1barcodeRC;
		protected boolean read2barcodeRC;
		
		/**
		 * @param mate1barcode Barcode that can be contained in read1
		 * @param mate2barcode Barcode that can be contained in read2
		 * @param maxMismatchMate1barcode Max mismatch in read1 barcode
		 * @param maxMismatchMate2barcode Max mismatch in read2 barcode
		 * @param mate1barcodeRC Whether read1 barcode appears reverse complemented in read1 sequence
		 * @param mate2barcodeRC Whether read2 barcode appears reverse complemented in read2 sequence
		 */
		public EitherMateContainsBarcode(String mate1barcode, String mate2barcode, int maxMismatchMate1barcode, int maxMismatchMate2barcode, boolean mate1barcodeRC, boolean mate2barcodeRC) {
			mateSpecificBarcode = new HashMap<Mate, String>();
			mateSpecificBarcode.put(Mate.MATE1, mate1barcode);
			mateSpecificBarcode.put(Mate.MATE2, mate2barcode);
			maxMismatchRead1barcode = maxMismatchMate1barcode;
			maxMismatchRead2barcode = maxMismatchMate2barcode;
			read1barcodeRC = mate1barcodeRC;
			read2barcodeRC = mate2barcodeRC;
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
		
		/**
		 * @return Barcode to look for in mate 1: the mate 1 barcode or its reverse complement depending on the boolean parameter
		 */
		protected String mate1barcodeOrRC() {
			if(read1barcodeRC) {
				return Sequence.reverseSequence(mateSpecificBarcode.get(Mate.MATE1));
			} else return mateSpecificBarcode.get(Mate.MATE1);
		}
		
		/**
		 * @return Barcode to look for in mate 2: the mate 2 barcode or its reverse complement depending on the boolean parameter
		 */
		protected String mate2barcodeOrRC() {
			if(read2barcodeRC) {
				return Sequence.reverseSequence(mateSpecificBarcode.get(Mate.MATE2));
			} else return mateSpecificBarcode.get(Mate.MATE2);
		}
		
		@Override
		public String identifyBarcode(FastqSequence read1, FastqSequence read2) {
			boolean read1hasSpecificBarcode = readHasBarcode(read1, mate1barcodeOrRC(), maxMismatchRead1barcode);
			boolean read2hasSpecificBarcode = readHasBarcode(read2, mate2barcodeOrRC(), maxMismatchRead2barcode);
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

		@Override
		public FastqSequence trimRead1(FastqSequence read1, String barcode) {
			throw new UnsupportedOperationException("Not implemented");
		}

		@Override
		public FastqSequence trimRead2(FastqSequence read2, String barcode) {
			throw new UnsupportedOperationException("Not implemented");
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
			return (read1hasBarcode && trimReads) ? trimRead1(read1, barcode) : read1;
		}

		@Override
		public FastqSequence getRead2SeqToWrite(FastqSequence read2, String barcode) {
			return (!read1hasBarcode && trimReads) ? trimRead2(read2, barcode) : read2;
		}

		@Override
		public FastqSequence trimRead1(FastqSequence read1, String barcode) {
			throw new UnsupportedOperationException("Not implemented");
		}

		@Override
		public FastqSequence trimRead2(FastqSequence read2, String barcode) {
			throw new UnsupportedOperationException("Not implemented");
		}

	}
	
	/**
	 * The options for barcode matcher implementations
	 * @author prussell
	 *
	 */
	private enum MatchImplementation {
		
		PERFECT_MATCH,
		EITHER_MATE_CONTAINS_BARCODE,
		ONE_MATE_CONTAINS_BARCODES_AT_KNOWN_POSITIONS;
		
		public String toString() {
			switch(this) {
			case EITHER_MATE_CONTAINS_BARCODE:
				return "either_mate_contains_barcode";
			case PERFECT_MATCH:
				return "perfect_match";
			case ONE_MATE_CONTAINS_BARCODES_AT_KNOWN_POSITIONS:
				return "one_mate_contains_barcodes_at_known_positions";
			default:
				throw new IllegalArgumentException("Not implemented");
			}
		}
		
		public static MatchImplementation fromString(String s) {
			if(s.equals(PERFECT_MATCH.toString())) {
				return PERFECT_MATCH;
			}
			if(s.equals(EITHER_MATE_CONTAINS_BARCODE.toString())) {
				return EITHER_MATE_CONTAINS_BARCODE;
			}
			if(s.equals(ONE_MATE_CONTAINS_BARCODES_AT_KNOWN_POSITIONS.toString())) {
				return ONE_MATE_CONTAINS_BARCODES_AT_KNOWN_POSITIONS;
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
	 * @throws IOException 
	 */
	private static BarcodeMatcher getBarcodeMatcher(CommandLineParser p) throws IOException {
		String perfectMatchBarcodeFile = p.getStringArg(perfectMatchBarcodeFileOption);
		String barcodeImpl = p.getStringArg(barcodeMatcherOption);
		boolean perfectMatchRead1hasBarcode = p.getBooleanArg(perfectMatchBarcodeOnRead1Option);
		String oneMateKnownPositionsBarcodeFile = p.getStringArg(oneMateKnownPositionsFileOption);
		String oneMateKnownPositionsMate = p.getStringArg(oneMateKnownPositionsMateOption);
		int oneMateKnownPositionsMaxMismatch = p.getIntArg(oneMateKnownPositionsMaxMismatchOption);
		
		MatchImplementation implementation = MatchImplementation.fromString(barcodeImpl);
		
		switch(implementation) {
		case EITHER_MATE_CONTAINS_BARCODE:
			throw new UnsupportedOperationException("Not implemented");
		case ONE_MATE_CONTAINS_BARCODES_AT_KNOWN_POSITIONS:
			if(oneMateKnownPositionsBarcodeFile == null) {
				throw new IllegalArgumentException("Must provide " + oneMateKnownPositionsFileOption);
			}
			if(oneMateKnownPositionsMate == null) {
				throw new IllegalArgumentException("Must provide " + oneMateKnownPositionsMateOption);
			}
			if(oneMateKnownPositionsMaxMismatch < 0) {
				throw new IllegalArgumentException("Must provide " + oneMateKnownPositionsMaxMismatchOption);
			}
			return new PairedEndBarcodeSplitter().new OneMateContainsBarcodesAtKnownPositions(Mate.fromString(oneMateKnownPositionsMate), oneMateKnownPositionsBarcodeFile, oneMateKnownPositionsMaxMismatch);
		case PERFECT_MATCH:
			if(perfectMatchBarcodeFile == null) {
				throw new IllegalArgumentException("Must provide " + perfectMatchBarcodeFileOption);
			}
			return new PairedEndBarcodeSplitter().new PerfectMatch(perfectMatchBarcodeFile, perfectMatchRead1hasBarcode);
		default:
			throw new IllegalArgumentException("Not implemented: " + implementation.toString());
		}
		
	}
	
	private static String perfectMatchBarcodeFileOption = "-pmb";
	private static String read1fqOption = "-f1";
	private static String read2fqOption = "-f2";
	private static String perfectMatchBarcodeOnRead1Option = "-pmb1";
	private static String trimReadsOption = "-rb";
	private static String barcodeMatcherOption = "-bmi";
	private static String eitherMateContainsBarcodeMate1BarcodeOption = "-emb1";
	private static String eitherMateContainsBarcodeMate2BarcodeOption = "-emb2";
	private static String eitherMateContainsBarcodeMaxMismatchRead1Option = "-emm1";
	private static String eitherMateContainsBarcodeMaxMismatchRead2Option = "-emm2";
	private static String oneMateKnownPositionsFileOption = "-omkpbf";
	private static String oneMateKnownPositionsMateOption = "-omkpm";
	private static String oneMateKnownPositionsMaxMismatchOption = "-omkpmm";
	
	/**
	 * Get the part of the read before an ungapped match of the barcode
	 * Provide a supposed match start obtained by a fast method (e.g. bitap, which allows gaps), 
	 * then this method verifies that there is in fact an ungapped match there
	 * @param read Read
	 * @param barcode Barcode
	 * @param matchStart Supposed start of barcode match
	 * @param maxMismatches Max mismatches
	 * @return The part of the read before the match
	 */
	@SuppressWarnings("unused")
	private static FastqSequence getSequenceBeforeUngappedMatch(FastqSequence read, String barcode, int matchStart, int maxMismatches) {
		if(!LocalAlignmentUtils.containsFullLengthUngappedMatchAtPosition(read.getSequence(), barcode, matchStart, maxMismatches, true)) {
			throw new IllegalArgumentException("No ungapped match of " + barcode + " to " + read.getSequence() + " at position " + 
					matchStart + " with at most " + maxMismatches + " mismatches.");
		}
		return read.trimEndBases(read.getLength() - matchStart);
	}
	
	/**
	 * Get the part of the read after an ungapped match of the barcode
	 * Provide a supposed match start obtained by a fast method (e.g. bitap, which allows gaps), 
	 * then this method verifies that there is in fact an ungapped match there
	 * @param read Read
	 * @param barcode Barcode
	 * @param matchStart Supposed start of barcode match
	 * @param maxMismatches Max mismatches
	 * @return The part of the read after the match
	 */
	private static FastqSequence getSequenceAfterUngappedMatch(FastqSequence read, String barcode, int matchStart, int maxMismatches) {
		if(!LocalAlignmentUtils.containsFullLengthUngappedMatchAtPosition(read.getSequence(), barcode, matchStart, maxMismatches, true)) {
			throw new IllegalArgumentException("No ungapped match of " + barcode + " to " + read.getSequence() + " at position " + 
					matchStart + " with at most " + maxMismatches + " mismatches.");
		}
		return read.trimStartBases(matchStart + barcode.length());
	}
	
	/**
	 * Get the part of the read between two ungapped matches of the barcode
	 * Provide two supposed match starts obtained by a fast method (e.g. bitap, which allows gaps), 
	 * then this method verifies that there are in fact ungapped matches there
	 * @param read Read
	 * @param barcode Barcode
	 * @param match1Start Supposed start of barcode match 1
	 * @param match2Start Supposed start of barcode match 2
	 * @param maxMismatches Max mismatches
	 * @return The part of the read between the barcode matches
	 */
	private static FastqSequence getSequenceBetweenTwoUngappedMatches(FastqSequence read, String barcode, int match1start, int match2start, int maxMismatches) {
		// Check that the two matches exist
		if(!LocalAlignmentUtils.containsFullLengthUngappedMatchAtPosition(read.getSequence(), barcode, match1start, maxMismatches, true)) {
			throw new IllegalArgumentException("No ungapped match of " + barcode + " to " + read.getSequence() + " at position " + 
					match1start + " with at most " + maxMismatches + " mismatches.");
		}
		if(!LocalAlignmentUtils.containsFullLengthUngappedMatchAtPosition(read.getSequence(), barcode, match2start, maxMismatches, true)) {
			throw new IllegalArgumentException("No ungapped match of " + barcode + " to " + read.getSequence() + " at position " + 
					match2start + " with at most " + maxMismatches + " mismatches.");
		}
		// Check that the matches don't overlap
		if(match1start + barcode.length() >= match2start) {
			throw new IllegalArgumentException("No sequence between the two matches of " + barcode + " to " + read.getSequence() + 
					" (" + match1start + " and " + match2start + ")");
		}
		// Take the part of the read between the two matches
		int numToTrimFromBeginning = match1start + barcode.length();
		int numToTrimFromEnd = read.getLength() - match2start;
		return read.trimStartBases(numToTrimFromBeginning).trimEndBases(numToTrimFromEnd);
	}
	
	/**
	 * Get the part of the read between two ungapped matches of the barcode, or after the only match if there is only one
	 * Provide supposed match starts obtained by a fast method (e.g. bitap, which allows gaps), 
	 * then this method verifies that there are in fact ungapped matches there
	 * Throws an exception if there are neither 1 nor 2 matches
	 * @param read Read
	 * @param barcode Barcode
	 * @param matchStarts Supposed match starts
	 * @param maxMismatches Max mismatches
	 * @return The part of the read between the two barcode matches or after the one
	 */
	@SuppressWarnings("unused")
	private static FastqSequence getSeqBetweenTwoUngappedMatchesOrAfterOneMatch(FastqSequence read, String barcode, List<Integer> matchStarts, int maxMismatches) {
		if(matchStarts.isEmpty()) {
			throw new IllegalStateException("No matches of " + barcode + " to read " + read.getSequence());
		}
		Iterator<Integer> iter = matchStarts.iterator();
		int match1 = iter.next().intValue();
		// Pair has been identified as containing barcode
		if(matchStarts.size() > 2) {
			throw new IllegalStateException("Too many matches (>2) of barcode " + barcode + " to read " + read.getSequence());
		}
		if(matchStarts.size() == 1) {
			// Take the part of the read after the first match of barcode
			return getSequenceAfterUngappedMatch(read, barcode, match1, maxMismatches);
		}
		// Take the part of the read between the two matches of barcode
		int match2 = iter.next().intValue();
		return getSequenceBetweenTwoUngappedMatches(read, barcode, match1, match2, maxMismatches);
	}
	
	
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
		p.addBooleanArg(trimReadsOption, "Trim reads according to barcode matching implementation", false, true);
		p.addStringArg(barcodeMatcherOption, "Barcode matcher implementation. Options: " + MatchImplementation.commaSeparatedList(), true);
		p.addStringArg(eitherMateContainsBarcodeMate1BarcodeOption, "For either mate implementation, read1 barcode", false, null);
		p.addStringArg(eitherMateContainsBarcodeMate2BarcodeOption, "For either mate implementation, read2 barcode", false, null);
		p.addIntArg(eitherMateContainsBarcodeMaxMismatchRead1Option, "For either mate implementation, max mismatches in read 1 barcode", false, -1);
		p.addIntArg(eitherMateContainsBarcodeMaxMismatchRead2Option, "For either mate implementation, max mismatches in read 2 barcode", false, -1);
		p.addStringArg(oneMateKnownPositionsFileOption, "For one mate known positions implementation, file of barcode and position in read (line format: barcode   position)", false, null);
		p.addStringArg(oneMateKnownPositionsMateOption, "For one mate known positions implementation, mate containing barcode. Options: " + Mate.MATE1.toString() + "," + Mate.MATE2.toString(), false, null);
		p.addIntArg(oneMateKnownPositionsMaxMismatchOption, "For one mate known positions implementation, max mismatches in barcode", false, -1);
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
