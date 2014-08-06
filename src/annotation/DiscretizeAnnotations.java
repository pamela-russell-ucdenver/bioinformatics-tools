/**
 * 
 */
package annotation;

import general.CommandLineParser;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import broad.pda.annotation.BEDFileParser;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Gene;
import nextgen.core.utils.AnnotationUtils;

/**
 * @author prussell
 * Discretize each annotation into regions that do not overlap another set; regions that overlap another set
 */
public class DiscretizeAnnotations {
	
	private Map<String, Collection<Gene>> regions;
	private Map<String, Collection<Gene>> otherRegions;
	private static Logger logger = Logger.getLogger(DiscretizeAnnotations.class.getName());
	
	private DiscretizeAnnotations(String bedFileRegions, String bedFileOtherRegions) throws IOException {
		this(BEDFileParser.loadDataByChr(bedFileRegions), BEDFileParser.loadDataByChr(bedFileOtherRegions));
	}
	
	private DiscretizeAnnotations(Map<String, Collection<Gene>> regionsByChr, Map<String, Collection<Gene>> otherRegionsByChr) {
		regions = regionsByChr;
		otherRegions = otherRegionsByChr;
	}
	
	private Annotation getNonOverlapAsAnnotation(Gene gene) {
		Collection<Gene> overlappers = AnnotationUtils.getChildren(gene, otherRegions);
		Gene rtrn = gene.copy();
		rtrn.setName(gene.getName() + "_minus_overlappers");
		if(overlappers.isEmpty()) {
			return rtrn;
		}
		return rtrn.minus(overlappers);
	}
	
	private Annotation getOverlapAsAnnotation(Gene gene) {
		Collection<Gene> overlappers = AnnotationUtils.getChildren(gene, otherRegions);
		if(overlappers == null) {
			logger.debug(gene.toUCSC() + " overlappers is null");
			return null;
		}
		if(overlappers.isEmpty()) {
			logger.debug(gene.toUCSC() + " overlappers is empty");
			return null;
		}
		Gene rtrn = new Gene(overlappers, gene.getName() + "_overlappers_only");
		return rtrn;
	}
	
	
	private void writeOverlap(String outBedFile) throws IOException {
		logger.info("Writing overlap to file " + outBedFile + "...");
		FileWriter w = new FileWriter(outBedFile);
		for(String chr : regions.keySet()) {
			logger.info(chr);
			if(!regions.containsKey(chr)) {
				continue;
			}
			for(Gene gene : regions.get(chr)) {
				Annotation overlap = getOverlapAsAnnotation(gene);
				if(overlap == null) {
					//logger.debug("Gene " + gene.getName() + " has no overlappers.");
					continue;
				}
				logger.debug("Overlap as annotation: " + overlap.toBED());
				w.write(overlap.toBED() + "\n");
			}
		}
		w.close();
		logger.info("Done writing file.");
	}
	
	private void writeNonOverlap(String outBedFile) throws IOException {
		logger.info("Writing non-overlap to file " + outBedFile + "...");
		FileWriter w = new FileWriter(outBedFile);
		for(String chr : regions.keySet()) {
			logger.info(chr);
			if(!regions.containsKey(chr)) {
				continue;
			}
			for(Gene gene : regions.get(chr)) {
				Annotation nonOverlap = getNonOverlapAsAnnotation(gene);
				if(nonOverlap == null) {
					//logger.debug("Gene " + gene.getName() + " has no overlappers.");
					continue;
				}
				logger.debug("Non-overlap as annotation: " + nonOverlap.toBED());
				w.write(nonOverlap.toBED() + "\n");
			}
		}
		w.close();
		logger.info("Done writing file.");
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-b", "Bed file of genes", true);
		p.addStringArg("-r", "Bed file of child regions", true);
		p.addStringArg("-o", "Output bed file prefix", true);
		p.addBooleanArg("-d", "Debug logging", false, false);
		p.parse(args);
		boolean debug = p.getBooleanArg("-d");
		if(debug) {
			logger.setLevel(Level.DEBUG);
			AnnotationUtils.logger.setLevel(Level.DEBUG);
		}
		String bedFile = p.getStringArg("-b");
		String regionFile = p.getStringArg("-r");
		String outPrefix = p.getStringArg("-o");
		
		DiscretizeAnnotations d = new DiscretizeAnnotations(bedFile, regionFile);
		d.writeOverlap(outPrefix + "_overlappers.bed");
		d.writeNonOverlap(outPrefix + "_non_overlappers.bed");
		
		logger.info("");
		logger.info("All done.");

	}

}
