package bam;

import guttmanlab.core.util.CommandLineParser;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;

import org.apache.log4j.Logger;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFormatException;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import broad.core.math.EmpiricalDistribution;

/**
 * @author prussell
 *
 */
public class ReadSizeDistribution {
	
	private String bamFile;
	private EmpiricalDistribution dist;
	private static Logger logger = Logger.getLogger(ReadSizeDistribution.class.getName());
	
	private ReadSizeDistribution(String bam) {
		bamFile = bam;
	}
	
	private void makeDistribution() {
		logger.info("Making distribution...");
		SAMFileReader reader = new SAMFileReader(new File(bamFile));
		SAMRecordIterator iter = reader.iterator();
		int numDone = 0;
		Collection<Double> readSizes = new ArrayList<Double>();
		
		while(iter.hasNext()) {
			try {
				SAMRecord record = iter.next();
				numDone++;
				if(numDone % 100000 == 0) {
					logger.info("Finished " + numDone + " records.");
				}
				readSizes.add(Double.valueOf(record.getReadLength()));
			} catch(SAMFormatException e) {
				logger.info("Skipping record: " + e.getMessage());
			}
		}
		
		reader.close();
		dist = new EmpiricalDistribution(readSizes);
		logger.info("Done making distribution.");
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-i", "Input sam or bam file", true);
		p.parse(args);
		String input = p.getStringArg("-i");
		
		ReadSizeDistribution r = new ReadSizeDistribution(input);
		r.makeDistribution();
		logger.info("Median read size is " + r.dist.getMedianOfAllDataValues());
		
		logger.info("All done.");

	}

}
