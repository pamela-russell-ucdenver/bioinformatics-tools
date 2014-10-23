/**
 * 
 */
package conservation;

import guttmanlab.core.util.CommandLineParser;
import guttmanlab.core.util.CountLogger;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeSet;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import broad.core.siphy.EvolutionaryModel.OmegaFit;
import broad.core.siphy.tools.conservation.EstimateOmegaPerExon;
import broad.pda.annotation.BEDFileParser;
import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Gene;

/**
 * @author prussell
 *
 */
public class Omega {
	
	private String alnDir;
	private File modFile;
	private static Logger logger = Logger.getLogger(Omega.class.getName());
	
	/**
	 * @param alignmentDir Directory containing maf files for each chromosome
	 * @param modelFile Evolutionary model file
	 */
	public Omega(String alignmentDir, String modelFile) {
		alnDir = alignmentDir;
		modFile = new File(modelFile);
	}
	
	/**
	 * Fit omega to a region
	 * @param region The region
	 * @return OmegaFit object
	 * @throws Exception
	 */
	public OmegaFit getOmegaFit(Annotation region) throws Exception {
		String alnFile = alnDir + "/" + region.getChr() + ".maf";
		try {
			return EstimateOmegaPerExon.fitOmega(region, modFile, alnFile, "maf");
		} catch (IllegalArgumentException e) {
			logger.warn("Caught exception when fitting omega for gene " + region.getName() + ". Returning null.");
			return null;
		}
	}
	
	/**
	 * Get omega score for a region
	 * @param region The region
	 * @return The omega score
	 * @throws Exception
	 */
	public double getOmega(Annotation region) throws Exception {
		OmegaFit of = getOmegaFit(region);
		if(of != null) {
			return of.getOmega();
		}
		return Double.NaN;
	}
	
	
	/**
	 * Write regions to bed file with fitted omega score in score field
	 * @param genesByChr Gene collection by chromosome
	 * @param outBed Output file
	 * @throws Exception
	 */
	public void writeOmegaAsBedScore(Map<String, Collection<Gene>> genesByChr, String outBed) throws Exception {
		writeOmegaAsBedScore(genesByChr, outBed, null);
	}
	
	/**
	 * Write regions to bed file with fitted omega score in score field
	 * @param genesByChr Gene collection by chromosome
	 * @param outBed Output file
	 * @param singleChr Single chromosome to write or null if all chromosomes
	 * @throws Exception
	 */
	public void writeOmegaAsBedScore(Map<String, Collection<Gene>> genesByChr, String outBed, String singleChr) throws Exception {
		FileWriter w = new FileWriter(outBed);
		logger.info("Writing omega scores as bed score to file " + outBed + "...");
		Collection<String> chrs = new TreeSet<String>();
		if(singleChr == null) {
			chrs.addAll(genesByChr.keySet());
		} else {
			chrs.add(singleChr);
		}
		int totalGenes = 0;
		for(String chr : chrs) {
			totalGenes += genesByChr.get(chr).size();
		}
		logger.info("Fitting omega for " + totalGenes + " regions...");
		CountLogger countLogger = new CountLogger(totalGenes, 100);
		for(String chr : chrs) {
			logger.info(chr);
			for(Gene gene : genesByChr.get(chr)) {
				countLogger.advance();
				gene.setBedScore(getOmega(gene));
				logger.debug(gene.getName() + "\t" + gene.getBedScore());
				w.write(gene.toBED() + "\n");
			}
		}
		w.close();
		logger.info("Done writing bed file.");
	}
	
	/**
	 * @param args
	 * @throws Exception 
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException, Exception {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-i", "Input bed file of regions", true);
		p.addStringArg("-a", "Alignment directory containing one maf file per chromosome", true);
		p.addStringArg("-m", "Evolutionary model file", true);
		p.addStringArg("-o", "Output bed file with fitted omega score in score field", false, null);
		p.addStringArg("-c", "Single chromosome to use", false, null);
		p.addBooleanArg("-d", "Debug logging", false, false);
		p.parse(args);
		boolean debug = p.getBooleanArg("-d");
		if(debug) {
			logger.setLevel(Level.DEBUG);
		}
		String inputBed = p.getStringArg("-i");
		String alignmentDir = p.getStringArg("-a");
		String modelFile = p.getStringArg("-m");
		String outputBed = p.getStringArg("-o");
		String singleChr = p.getStringArg("-c");
		
		Omega omega = new Omega(alignmentDir, modelFile);
		
		if(outputBed != null) {
			omega.writeOmegaAsBedScore(BEDFileParser.loadDataByChr(inputBed), outputBed, singleChr);
		}
	}

}
