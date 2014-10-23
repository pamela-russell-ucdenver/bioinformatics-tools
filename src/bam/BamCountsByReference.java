package bam;

import guttmanlab.core.util.CommandLineParser;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;

import org.apache.log4j.Logger;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFormatException;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class BamCountsByReference {
	
	private static Logger logger = Logger.getLogger(BamCountsByReference.class.getName());
	private String bam;
	private Map<String, Integer> countsByReference;

	private BamCountsByReference(String bamFile) {
		bam = bamFile;
		countsByReference = new HashMap<String, Integer>();
	}
	
	private void incrementCount(String refName) {
		if(!countsByReference.containsKey(refName))	{
			countsByReference.put(refName, Integer.valueOf(1));
			return;
		}
		int prevCount = countsByReference.get(refName).intValue();
		countsByReference.put(refName, Integer.valueOf(prevCount + 1));
	}
	
	private void writeToFile(String outFile) throws IOException {
		TreeMap<String, Integer> sortedMap = new TreeMap<String, Integer>();
		sortedMap.putAll(countsByReference);
		FileWriter w = new FileWriter(outFile);
		for(String ref : sortedMap.keySet()) {
			w.write(ref + "\t" + sortedMap.get(ref).toString() + "\n");
		}
		w.close();
	}
	
	private void makeCounts() {
		logger.info("Making counts...");
		SAMFileReader reader = new SAMFileReader(new File(bam));
		SAMRecordIterator iter = reader.iterator();
		int numDone = 0;
		
		while(iter.hasNext()) {
			try {
				SAMRecord record = iter.next();
				numDone++;
				if(numDone % 1000000 == 0) {
					logger.info("Finished " + numDone + " records.");
				}
				if(!record.getReadUnmappedFlag()) {
					incrementCount(record.getReferenceName());
				}
			} catch(SAMFormatException e) {
				logger.info("Skipping record: " + e.getMessage());
			}
		}
		
		reader.close();
		logger.info("Done making counts.");
	}
	
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-b", "Bam file", true);
		p.addStringArg("-o", "Output table", true);
		p.parse(args);
		String bam = p.getStringArg("-b");
		String out = p.getStringArg("-o");
		
		BamCountsByReference b = new BamCountsByReference(bam);
		b.makeCounts();
		b.writeToFile(out);
		
	}

}
