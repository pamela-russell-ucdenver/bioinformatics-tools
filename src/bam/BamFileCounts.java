/**
 * 
 */
package bam;

import general.CommandLineParser;

import java.io.File;

import org.apache.log4j.Logger;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFormatException;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;



/**
 * @author prussell
 *
 */
public class BamFileCounts {
	
	private String bamFile;
	private int mapped;
	private int unmapped;
	private static Logger logger = Logger.getLogger(BamFileCounts.class.getName());
	
	private BamFileCounts(String bam) {
		bamFile = bam;
		mapped = 0;
		unmapped = 0;
	}
	
	private void makeCounts() {
		logger.info("Making counts...");
		SAMFileReader reader = new SAMFileReader(new File(bamFile));
		SAMRecordIterator iter = reader.iterator();
		int numDone = 0;
		mapped = 0;
		unmapped = 0;
		
		while(iter.hasNext()) {
			try {
				SAMRecord record = iter.next();
				numDone++;
				if(numDone % 100000 == 0) {
					logger.info("Finished " + numDone + " records.");
				}
				if(!record.getReadUnmappedFlag()) mapped++;
				else unmapped++;
			} catch(SAMFormatException e) {
				logger.info("Skipping record: " + e.getMessage());
			}
		}
		
		reader.close();
		logger.info("Done making counts.");
	}


	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-b", "Sam or bam file", true);
		p.parse(args);
		String bam = p.getStringArg("-b");
		BamFileCounts c = new BamFileCounts(bam);
		c.makeCounts();
		int mapped = c.mapped;
		int unmapped = c.unmapped;
		int total = mapped + unmapped;
		double mappedPct = 100 * (double) mapped / total;
		double unmappedPct = 100 * (double) unmapped / total;
		logger.info("Total records:\t" + total);
		logger.info("Mapped:\t" + mapped + " (" + mappedPct + "%)");
		logger.info("Unmapped:\t" + unmapped + " (" + unmappedPct + "%)");

	}

}
