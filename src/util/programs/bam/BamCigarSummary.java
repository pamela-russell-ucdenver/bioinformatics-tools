package util.programs.bam;

import guttmanlab.core.util.CommandLineParser;

import java.io.File;
import java.util.List;

import org.apache.log4j.Logger;

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.TextCigarCodec;

public class BamCigarSummary {
	
	private long totalRecords;
	private long totalMappedRecords;
	private long totalUnmappedRecords;
	private long totalNumInsertions;
	private long totalInsertionLength;
	private long totalNumDeletions;
	private long totalDeletionLength;
	private long totalMappedLength;
	private SAMFileReader reader;
	private static Logger logger = Logger.getLogger(BamCigarSummary.class.getName());
	
	private BamCigarSummary(String bamFile) {
		reader = new SAMFileReader(new File(bamFile));
		totalRecords = 0;
		totalMappedRecords = 0;
		totalUnmappedRecords = 0;
		totalNumInsertions = 0;
		totalInsertionLength = 0;
		totalNumDeletions = 0;
		totalDeletionLength = 0;
		totalMappedLength = 0;
		processFile();
	}
	
	private void countRecord(SAMRecord record) {
		
		totalRecords++;
		
		if(record.getReadUnmappedFlag()) {
			totalUnmappedRecords++;
			return;
		}
		
		totalMappedRecords++;
		
    	Cigar cigar = TextCigarCodec.getSingleton().decode(record.getCigarString());
    	List<CigarElement> elements=cigar.getCigarElements();
		
		for(CigarElement element : elements){
			CigarOperator op = element.getOperator();
			int length = element.getLength();
			if(op.equals(CigarOperator.MATCH_OR_MISMATCH)){
				totalMappedLength += length;
			} else if(op.equals(CigarOperator.INSERTION)) {
				totalNumInsertions++;
				totalInsertionLength += length;
			} else if(op.equals(CigarOperator.DELETION)) {
				totalNumDeletions++;
				totalDeletionLength += length;
			}
		}

	}
	
	private void processFile() {
		logger.info("");
		logger.info("Processing bam file...");
		SAMRecordIterator iter = reader.iterator();
		long numDone = 0;
		while(iter.hasNext()) {
			countRecord(iter.next());
			numDone++;
			if(numDone % 10000000 == 0) {
				logger.info("Counted " + numDone + " records.");
			}
		}
		logger.info("Done processing file.");
	}
	
	private void printSummary() {
		
		double avgMappedLength = (double) totalMappedLength / (double) totalMappedRecords;
		double avgNumInsertions = (double) totalNumInsertions / (double) totalMappedRecords;
		double avgInsertionLength = (double) totalInsertionLength / (double) totalNumInsertions;
		double avgNumDeletions = (double) totalNumDeletions / (double) totalMappedRecords;
		double avgDeletionLength = (double) totalDeletionLength / (double) totalNumDeletions;
		
		logger.info("");
		logger.info("Total records:\t" + totalRecords);
		logger.info("Mapped records:\t" + totalMappedRecords);
		logger.info("Unmapped records:\t" + totalUnmappedRecords);
		logger.info("Average mapped length per mapped record:\t" + avgMappedLength);
		logger.info("Average insertions per mapped record:\t" + avgNumInsertions);
		logger.info("Average insertion length:\t" + avgInsertionLength);
		logger.info("Average deletions per mapped record:\t" + avgNumDeletions);
		logger.info("Average deletion length:\t" + avgDeletionLength);
	}
	
	public static void main(String[] args) {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-b", "Bam file", true);
		p.parse(args);
		String bamFile = p.getStringArg("-b");
		BamCigarSummary bcs = new BamCigarSummary(bamFile);
		bcs.printSummary();

	}

}
