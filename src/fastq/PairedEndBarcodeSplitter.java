package fastq;

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
	
	private String read1fq;
	private String read2fq;
	private boolean read1hasBarcode;
	private Map<String, String> read1out;
	private Map<String, String> read2out;
	private Map<String, BufferedWriter> read1writers;
	private Map<String, BufferedWriter> read2writers;
	private BufferedWriter unmatched1writer;
	private BufferedWriter unmatched2writer;
	private Collection<String> barcodes;
	private boolean removeBarcode;
	private FastqParser fastq1parser;
	private FastqParser fastq2parser;
	
	private static Logger logger = Logger.getLogger(PairedEndBarcodeSplitter.class.getName());

	/**
	 * @param read1fastq Read1 fastq file
	 * @param read2fastq Read2 fastq file
	 * @param barcodeFile File containing simple list of barcodes
	 * @param read1barcode True iff the barcode is in read1
	 * @param removeBarcodeFromReads Remove the barcode before writing
	 * @throws IOException 
	 */
	private PairedEndBarcodeSplitter(String read1fastq, String read2fastq, String barcodeFile, boolean read1barcode, boolean removeBarcodeFromReads) throws IOException {
		read1fq = read1fastq;
		read2fq = read2fastq;
		read1hasBarcode = read1barcode;
		barcodes = new HashSet<String>();
		barcodes.addAll(FileUtil.fileLinesAsList(barcodeFile));
		removeBarcode = removeBarcodeFromReads;
		resetReaders();
		resetWriters();
	}
	
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
	 * @param predicate Predicate for whether a read has the barcode
	 * @throws IOException
	 */
	private void writeMatches(ReadHasBarcode predicate) throws IOException {
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
			String barcode = identifyBarcode(read1, read2, predicate);
			if(barcode != null) {
				FastqSequence toWrite1 = (read1hasBarcode && removeBarcode) ? removeBarcode(read1, barcode) : read1;
				FastqSequence toWrite2 = (!read1hasBarcode && removeBarcode) ? removeBarcode(read2, barcode) : read2;
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
		// TODO
		throw new UnsupportedOperationException("not implemented");
	}
	
	/**
	 * @param read1 Read 1
	 * @param read2 Read 2
	 * @param predicate A predicate that is true iff a read has a match
	 * @return The barcode for the read pair or null if no match
	 */
	private String identifyBarcode(FastqSequence read1, FastqSequence read2, ReadHasBarcode predicate) {
		FastqSequence barcodeRead = read1hasBarcode ? read1 : read2;
		for(String barcode : barcodes) {
			if(predicate.readHasBarcode(barcodeRead, barcode)) {
				return barcode;
			}
		}
		return null;
	}
	
	/**
	 * Predicate representing whether a read contains a barcode
	 * @author prussell
	 *
	 */
	private interface ReadHasBarcode {
		public boolean readHasBarcode(FastqSequence read, String barcode);
	}
	
	/**
	 * Check whether a fastq record contains a perfect match of a sequence
	 * @author prussell
	 *
	 */
	private class ReadContainsPerfectMatch implements ReadHasBarcode {
		
		public ReadContainsPerfectMatch() {}
		
		@Override
		public boolean readHasBarcode(FastqSequence read, String barcode) {
			String readSeq = read.getSequence();
			return readSeq.contains(barcode);
		}
		
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-b", "File containing simple list of barcodes", true);
		p.addStringArg("-f1", "Read1 fastq file", true);
		p.addStringArg("-f2", "Read2 fastq file", true);
		p.addBooleanArg("-b1", "true if barcode is on read1, false if read2", true);
		p.addBooleanArg("-rb", "Remove barcode when writing separated files", false, true);
		p.parse(args);
		String read1fastq = p.getStringArg("-f1");
		String read2fastq = p.getStringArg("-f2");
		String barcodeFile = p.getStringArg("-b");
		boolean read1barcode = p.getBooleanArg("-b1");
		boolean removeBarcodeFromReads = p.getBooleanArg("-rb");
		
		PairedEndBarcodeSplitter bs = new PairedEndBarcodeSplitter(read1fastq, read2fastq, barcodeFile, read1barcode, removeBarcodeFromReads);
		
		bs.writeMatches(bs.new ReadContainsPerfectMatch());
		
		logger.info("");
		logger.info("All done.");

	}

}
