package bam;

import guttmanlab.core.util.CommandLineParser;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Collection;
import java.util.HashSet;

import org.apache.log4j.Logger;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class BamFilterByCustomTag {
	
	private String tag;
	private Collection<String> validTagValues;
	private static Logger logger = Logger.getLogger(BamFilterByCustomTag.class.getName());
	
	/**
	 * @param customTag The custom tag field to use (should start with "X")
	 * @param validTagListFile File containing one valid tag value per line
	 * @throws IOException
	 */
	private BamFilterByCustomTag(String customTag, String validTagListFile) throws IOException {
		tag = customTag;
		readValidTagValues(validTagListFile);
	}
	
	/**
	 * Check whether the record has a valid value for the SAM tag
	 * @param record Record
	 * @return True iff the record's value for the tag is in the collection of valid values
	 */
	private boolean hasValidTagValue(SAMRecord record) {
		return validTagValues.contains(record.getAttribute(tag));
	}
	
	/**
	 * Initialize the set of valid SAM tag values
	 * @param listFile File containing one valid value per line
	 * @throws IOException
	 */
	private void readValidTagValues(String listFile) throws IOException {
		validTagValues = new HashSet<String>();
		BufferedReader r = new BufferedReader(new FileReader(listFile));
		while(r.ready()) {
			validTagValues.add(r.readLine());
		}
		r.close();
	}
	
	/**
	 * Filter the bam file, keeping records with a valid value for the SAM tag
	 * @param inputBam Input bam file to filter
	 * @param outputBam Output filtered bam file
	 */
	private void filterBam(String inputBam, String outputBam) {
		logger.info("");
		logger.info("Writing filtered bam file...");
		SAMFileReader reader = new SAMFileReader(new File(inputBam));
		SAMFileWriter writer = new SAMFileWriterFactory().makeBAMWriter(reader.getFileHeader(), false, new File(outputBam));
		SAMRecordIterator iter = reader.iterator();
		int numDone = 0;
		int kept = 0;
		while(iter.hasNext()) {
			SAMRecord record = iter.next();
			numDone++;
			if(hasValidTagValue(record)) {
				kept++;
				writer.addAlignment(record);
			}
			if(numDone % 1000000 == 0) {
				logger.info("Finished " + numDone + " records. Kept " + kept);
			}
		}
		iter.close();
		reader.close();
		writer.close();
	}
	
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-t", "SAM tag e.g. XC", true);
		p.addStringArg("-v", "File containing list of valid tag values, one per line", true);
		p.addStringArg("-i", "Input bam file to filter for valid tag values", true);
		p.addStringArg("-o", "Output filtered bam file to write", true);
		p.parse(args);
		String tag = p.getStringArg("-t");
		String list = p.getStringArg("-v");
		String input = p.getStringArg("-i");
		String output = p.getStringArg("-o");
		
		BamFilterByCustomTag b = new BamFilterByCustomTag(tag, list);
		b.filterBam(input, output);
		
		logger.info("");
		logger.info("All done.");
		
	}
	
}
