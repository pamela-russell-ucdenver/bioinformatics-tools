package util.programs.bam;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import org.apache.log4j.Logger;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import fastq.FastqRecord;
import guttmanlab.core.util.CommandLineParser;

public class RevertBam {
	
	private static Logger logger = Logger.getLogger(RevertBam.class.getName());
	
	/**
	 * Enum for different ReadNameGenerator implementations
	 * @author prussell
	 *
	 */
	private static enum ReadNameGeneratorType {
		
		
		MCCARROLL_SINGLE_CELL, // McCarroll lab 2015 Cell paper Drop-seq
		QNAME; // Just the read name
		
		public String toString() {
			switch(this) {
			case MCCARROLL_SINGLE_CELL:
				return "McCarroll_single_cell";
			case QNAME:
				return ("Qname");
			default:
				throw new IllegalArgumentException("Not implemented");
			}
		}
		
		public static ReadNameGeneratorType fromString(String name) {
			if(name.equals(MCCARROLL_SINGLE_CELL.toString())) {
				return MCCARROLL_SINGLE_CELL;
			}
			if(name.equals(QNAME.toString())) {
				return QNAME;
			}
			throw new IllegalArgumentException("Not implemented");
		}
	
		public static String getCommaSeparatedList() {
			ReadNameGeneratorType[] vals = values();
			String rtrn = vals[0].toString();
			for(int i = 1; i < vals.length; i++) {
				rtrn += ", " + vals[i].toString();
			}
			return rtrn;
		}
		
	}
	
	/**
	 * Get ReadNameGeneratorType from string name
	 * @param readNameGeneratorType Name of implmentation
	 * @return The ReadNameGeneratorType represented by the name
	 */
	private static ReadNameGenerator readNameGeneratorFactory(String readNameGeneratorType) {
		ReadNameGeneratorType nameGeneratorType = ReadNameGeneratorType.fromString(readNameGeneratorType);
		switch(nameGeneratorType) {
		case MCCARROLL_SINGLE_CELL:
			return new McCarrollSingleCellReadName();
		case QNAME:
			return new Qname();
		default:
			throw new IllegalArgumentException("Not implemented");
		}
	}
	
	/**
	 * A way to generate a new name for a read, beyond just QNAME, from a SAMRecord
	 * @author prussell
	 *
	 */
	private interface ReadNameGenerator {
		
		/**
		 * Get a name for a SAM record
		 * @param record SAM record
		 * @return The name
		 */
		public String createName(SAMRecord record);
		
	}
	
	/**
	 * Just get the read name
	 * @author prussell
	 *
	 */
	private static final class Qname implements ReadNameGenerator {
		
		public Qname() {}
		
		@Override
		public String createName(SAMRecord record) {
			return record.getReadName();
		}
		
	}
	
	/**
	 * Format for McCarroll lab 2015 Cell paper Drop-seq
	 * @author prussell
	 *
	 */
	private static final class McCarrollSingleCellReadName implements ReadNameGenerator {
		
		public McCarrollSingleCellReadName() {}
		
		@Override
		public String createName(SAMRecord record) {
			String cellBarcode = record.getAttribute("XC").toString();
			String umi = record.getAttribute("XM").toString();
			return record.getReadName() + " " + cellBarcode + " " + umi;
		}
		
	}
	
	/**
	 * 
	 * @param record
	 * @param name
	 * @param writer
	 * @throws IOException
	 */
	private static void writeReverted(SAMRecord record, String name, BufferedWriter writer) throws IOException {
		FastqRecord fastqRecord = new FastqRecord(record, name);
		fastqRecord.write(writer);
	}
	
	
	private static void writeRevertedFastq(String bamFile, String outFastq, String readNameGeneratorType) throws IOException {
		logger.info("");
		logger.info("Reverting " + bamFile + " to fastq format and writing to " + outFastq + "...");
		int numDone = 0;
		BufferedWriter writer = new BufferedWriter(new FileWriter(outFastq));
		ReadNameGenerator nameGenerator = readNameGeneratorFactory(readNameGeneratorType);
		SAMFileReader reader = new SAMFileReader(new File(bamFile));
		SAMRecordIterator iter = reader.iterator();
		while(iter.hasNext()) {
			SAMRecord samRecord = iter.next();
			writeReverted(samRecord, nameGenerator.createName(samRecord), writer);
			numDone++;
			if(numDone % 1000000 == 0) {
				logger.info("Finished " + numDone + " records.");
			}
		}
		iter.close();
		writer.close();
	}
	
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-i", "Input bam", true);
		p.addStringArg("-o", "Output fastq", true);
		p.addStringArg("-n", "Read name generation scheme (options = " + ReadNameGeneratorType.getCommaSeparatedList() + ")", true);
		p.parse(args);
		String input = p.getStringArg("-i");
		String output = p.getStringArg("-o");
		String nameGenerator = p.getStringArg("-n");
		
		writeRevertedFastq(input, output, nameGenerator);
		
		logger.info("");
		logger.info("All done.");
		
	}

}
