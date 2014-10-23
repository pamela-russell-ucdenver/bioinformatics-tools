package bam;

import guttmanlab.core.util.CommandLineParser;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import org.apache.log4j.Logger;

import net.sf.samtools.BAMFileWriter;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFormatException;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMValidationError;

/**
 * @author prussell
 *
 */
public class ValidateBam {
	
	private static Logger logger = Logger.getLogger(ValidateBam.class.getName());
	
	private static void validateFile(String input, String output) {
		
		logger.info("Reading file " + input + " and writing valid records to file " + output + "...");
		
		SAMFileReader reader = new SAMFileReader(new File(input));
		reader.setValidationStringency(SAMFileReader.ValidationStringency.STRICT);
		SAMFileHeader header = reader.getFileHeader();
		SAMRecordIterator iter = reader.iterator();
		
		BAMFileWriter writer = new BAMFileWriter(new File(output));
		writer.setHeader(header);
		
		int numDone = 0;
		
		while(iter.hasNext()) {
			try {
				SAMRecord record = iter.next();
				numDone++;
				if(numDone % 100000 == 0) {
					logger.info("Finished " + numDone + " records.");
				}
				
				// Check that cigar length matches read length
				Cigar cigar = record.getCigar();
				int cigarLength = 0;
				for(CigarElement e : cigar.getCigarElements()) {
					CigarOperator o = e.getOperator();
					if(o.equals(CigarOperator.M) || o.equals(CigarOperator.I) || o.equals(CigarOperator.S) || o.equals(CigarOperator.EQ) || o.equals(CigarOperator.X)) {
						cigarLength += e.getLength();
					}
				}
				int readLength = cigar.getReadLength();
				if(cigarLength != readLength) {
					logger.error("Skipping record " + record.getReadName() + " because cigar length (" + cigarLength + ") does not equal read length (" + readLength + ").");
					continue;
				}
				
				// Check other errors
				List<SAMValidationError> errors = record.isValid();
				List<SAMValidationError> cigarErrors = record.validateCigar(-1);
				List<SAMValidationError> allErrors = new ArrayList<SAMValidationError>();
				if(errors != null) {
					allErrors.addAll(errors);
				}
				if(cigarErrors != null) {
					allErrors.addAll(cigarErrors);
				}
				if(!allErrors.isEmpty()) {
					logger.error("Skipping record " + record.getReadName());
					for(SAMValidationError error : allErrors) {
						logger.error(error.getType() + "\t" + error.getMessage());
					}
					continue;
				}
				writer.addAlignment(record);
			} catch(SAMFormatException e) {
				logger.info("Skipping record: " + e.getMessage());
			}
		}
		
		reader.close();
		writer.close();
		
		logger.info("Done writing file.");
		
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-i", "Input sam or bam file", true);
		p.addStringArg("-o", "Output bam file", true);
		p.parse(args);
		String input = p.getStringArg("-i");
		String output = p.getStringArg("-o");
		
		validateFile(input, output);
		logger.info("All done.");
		
	}
	
}
