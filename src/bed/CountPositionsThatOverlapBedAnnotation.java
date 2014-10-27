package bed;

import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.FeatureCollection;
import guttmanlab.core.util.CommandLineParser;
import guttmanlab.core.util.StringParser;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Map;

import org.apache.log4j.Logger;

public class CountPositionsThatOverlapBedAnnotation {
	
	private Map<String, FeatureCollection<Gene>> genesByChr;
	private static Logger logger = Logger.getLogger(CountPositionsThatOverlapBedAnnotation.class.getName());
	private int overlappers;
	private int nonoverlappers;
	
	/**
	 * @param bedFile Gene annotation in bed format
	 * @param referenceSizes Reference genome size file. Line format: chr   size
	 * @param positionList File of positions to check for gene overlap. Line format: chr   pos
	 * @throws IOException
	 */
	private CountPositionsThatOverlapBedAnnotation(String bedFile, String referenceSizes, String positionList) throws IOException {
		loadGenes(bedFile, referenceSizes);
		countOverlappers(positionList);
	}
	
	private void loadGenes(String bedFile, String referenceSizes) throws IOException {
		genesByChr = BEDFileIO.loadFromFileByReferenceName(bedFile, referenceSizes);
	}
	
	/**
	 * Count number of positions in list that overlap some exon from annotation
	 * @param positionList List file. Line format: chr   pos
	 * @throws IOException 
	 */
	private void countOverlappers(String positionList) throws IOException {
		logger.info("Counting overlappers...");
		overlappers = 0;
		nonoverlappers = 0;
		BufferedReader reader = new BufferedReader(new FileReader(positionList));
		StringParser p = new StringParser();
		while(reader.ready()) {
			p.parse(reader.readLine());
			int numFields = p.getFieldCount();
			if(numFields == 0) continue;
			if(numFields != 2) {
				reader.close();
				throw new IllegalArgumentException("Line format for position file: chr   pos");
			}
			String chr = p.asString(0);
			int pos = p.asInt(1);
			SingleInterval interval = new SingleInterval(chr, pos, pos + 1);
			if(!genesByChr.containsKey(chr)) {
				//logger.warn("Skipping position " + chr + " " + pos + " because there are no annotations on the chromosome.");
				continue;
			}
			if(genesByChr.get(chr).overlaps(interval)) {
				overlappers++;
				logger.info(chr + "\t" + pos);
			} else {
				nonoverlappers++;
			}
		}
		reader.close();
		logger.info("Done counting overlappers.");
	}
	
	private int getNumOverlappers() {
		return overlappers;
	}
	
	private int getNumNonOverlappers() {
		return nonoverlappers;
	}
	
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-b", "Bed gene annotation", true);
		p.addStringArg("-c", "Chromosome size file", true);
		p.addStringArg("-p", "Position list file", true);
		p.parse(args);
		String bedFile = p.getStringArg("-b");
		String referenceSizes = p.getStringArg("-c");
		String positionList = p.getStringArg("-p");
		
		CountPositionsThatOverlapBedAnnotation c = new CountPositionsThatOverlapBedAnnotation(bedFile, referenceSizes, positionList);
		
		int overlappers = c.getNumOverlappers();
		int nonoverlappers = c.getNumNonOverlappers();
		
		logger.info("There are " + overlappers + " overlappers and " + nonoverlappers + " non-overlappers.");
		logger.info("");
		logger.info("All done.");
		
	}

}
