package bam;

import guttmanlab.core.util.CommandLineParser;
import guttmanlab.core.util.StringParser;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
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

public class AlignmentCountsByReference {
	
	private static Logger logger = Logger.getLogger(AlignmentCountsByReference.class.getName());
	private String alignmentFile;
	private Map<String, Integer> countsByReference;

	private AlignmentCountsByReference(String input) {
		alignmentFile = input;
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
	
	private void makeCountsMap() throws IOException {
		logger.info("Making counts...");
		BufferedReader r = new BufferedReader(new FileReader(alignmentFile));
		int numDone = 0;
		
		StringParser s = new StringParser();
		
		while(r.ready()) {
			String line = r.readLine();
			numDone++;
			if(numDone % 1000000 == 0) {
				logger.info("Finished " + numDone + " records.");
			}
			s.parse(line);
			incrementCount(s.asString(2));
		}
		
		r.close();
		logger.info("Done making counts.");
	}
	
	private void makeCountsBam(boolean primaryAlignmentsOnly) {
		logger.info("Making counts...");
		SAMFileReader reader = new SAMFileReader(new File(alignmentFile));
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
					if(primaryAlignmentsOnly && record.getNotPrimaryAlignmentFlag()) {
						continue;
					}
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
		p.addStringArg("-b", "Bam file", false, null);
		p.addStringArg("-m", "Map alignment file", false, null);
		p.addStringArg("-o", "Output table", true);
		p.addBooleanArg("-p", "For bam file, count primary alignments only", false, false);
		p.parse(args);
		String bam = p.getStringArg("-b");
		String out = p.getStringArg("-o");
		String map = p.getStringArg("-m");
		boolean primaryOnly = p.getBooleanArg("-p");
		
		if((bam == null && map == null) || (bam != null && map != null)) {
			throw new IllegalArgumentException("Provide one: -b or -m");
		}
		
		String alignmentFile = bam != null ? bam : map;
		
		AlignmentCountsByReference b = new AlignmentCountsByReference(alignmentFile);
		
		if(bam != null) {
			b.makeCountsBam(primaryOnly);
		} else {
			b.makeCountsMap();
		}
		
		b.writeToFile(out);
		
	}

}
