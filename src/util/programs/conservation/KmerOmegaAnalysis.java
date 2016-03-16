/**
 * 
 */
package util.programs.conservation;

import guttmanlab.core.util.CommandLineParser;
import guttmanlab.core.util.CountLogger;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeSet;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import util.OmegaKmerFileReader;
import broad.core.siphy.EvolutionaryModel.OmegaFit;
import broad.pda.annotation.BEDFileParser;
import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Annotation.Strand;
import nextgen.core.annotation.Gene;
import nextgen.core.utils.AnnotationUtils;

/**
 * @author prussell
 *
 */
public class KmerOmegaAnalysis {
	
	private OmegaKmerFileReader kmerFileReader;
	private Map<String, Collection<Gene>> genes;
	private static Logger logger = Logger.getLogger(KmerOmegaAnalysis.class.getName());
	
	/**
	 * @param omegaDirectory
	 * @param kmerSize
	 * @param fileChunkSize
	 * @param geneBedFile
	 * @param useCache
	 * @throws IOException
	 */
	public KmerOmegaAnalysis(String omegaDirectory, int kmerSize, int fileChunkSize, String geneBedFile, boolean useCache) throws IOException {
		this(new OmegaKmerFileReader(omegaDirectory, kmerSize, fileChunkSize, useCache), BEDFileParser.loadDataByChr(geneBedFile));
	}
	
	/**
	 * @param k
	 * @param genesByChr
	 */
	public KmerOmegaAnalysis(OmegaKmerFileReader k, Map<String, Collection<Gene>> genesByChr) {
		kmerFileReader = k;
		genes = genesByChr;
	}
	
	/**
	 * @param region The region
	 * @param pvalCutoff P value cutoff
	 * @return Merged significant kmers
	 * @throws IOException
	 */
	public Collection<Annotation> getMergedSignificantKmers(Gene region, double pvalCutoff) throws IOException {
		Map<Integer, OmegaFit> scores = kmerFileReader.getAllKmerScores(region);
		TreeSet<Annotation> sigKmers = new TreeSet<Annotation>();
		for(OmegaFit of : scores.values()) {
			if(of.getPVal() < pvalCutoff) {
				Annotation kmer = of.getRegion();
				kmer.setOrientation(Strand.POSITIVE);
				sigKmers.add(kmer);
			}
		}
		return AnnotationUtils.mergeOverlappingBlocks(sigKmers);
	}
	
	@SuppressWarnings("unused")
	private void writeOmegaQuantiles(double quantile, String outFile) throws IOException {
		writeOmegaQuantiles(quantile, outFile, null);
	}
		
	private void writeOmegaQuantiles(double quantile, String outFile, String singleChr) throws IOException {
		logger.info("Writing " + quantile + " omega quantile for each gene to file " + outFile + "...");
		FileWriter w = new FileWriter(outFile);
		int totalGenes = 0;
		Collection<String> chrsToUse = new TreeSet<String>();
		if(singleChr != null) {
			chrsToUse.add(singleChr);
		} else {
			chrsToUse.addAll(genes.keySet());
		}
		for(String chr : chrsToUse) {
			totalGenes += genes.get(chr).size();
		}
		CountLogger countLogger = new CountLogger(totalGenes, 100);
		for(String chr : chrsToUse) {
			logger.info(chr);
			for(Gene gene : genes.get(chr)) {
				double omegaQuantile = kmerFileReader.getOmegaQuantile(gene, quantile);
				w.write(gene.getName() + "\t" + omegaQuantile + "\n");
				countLogger.advance();
			}
		}
		w.close();
		logger.info("Done writing file.");
	}
	
	
	@SuppressWarnings("unused")
	private void writeMergedSignificantKmers(double pvalCutoff, String outBed) throws IOException {
		writeMergedSignificantKmers(pvalCutoff, outBed, null);
	}
	
	private void writeMergedSignificantKmers(double pvalCutoff, String outBed, String singleChr) throws IOException {
		logger.info("Writing merged kmers with P value < " + pvalCutoff + " to file " + outBed + "...");
		FileWriter w = new FileWriter(outBed);
		int totalGenes = 0;
		Collection<String> chrsToUse = new TreeSet<String>();
		if(singleChr != null) {
			chrsToUse.add(singleChr);
		} else {
			chrsToUse.addAll(genes.keySet());
		}
		for(String chr : chrsToUse) {
			totalGenes += genes.get(chr).size();
		}
		CountLogger countLogger = new CountLogger(totalGenes, 100);
		for(String chr : chrsToUse) {
			logger.info(chr);
			for(Gene gene : genes.get(chr)) {
				Collection<Annotation> sigRegions = getMergedSignificantKmers(gene, pvalCutoff);
				for(Annotation region : sigRegions) {
					region.setName(gene.getName() + ":" + region.toUCSC());
					w.write(region.toBED() + "\n");
				}
				countLogger.advance();
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
		p.addStringArg("-w", "Omega directory", true);
		p.addIntArg("-k", "Kmer size", true);
		p.addIntArg("-c", "Omega file chunk size", true);
		p.addStringArg("-g", "Gene bed file", true);
		p.addBooleanArg("-d", "Debug logging", false, false);
		p.addDoubleArg("-q", "Omega quantile", false, 0.9);
		p.addStringArg("-oq", "Output table for omega quantile of each gene", false, null);
		p.addStringArg("-os", "Output bed file of merged significant kmers", false, null);
		p.addDoubleArg("-p", "P value cutoff for significant kmer", false, 0.01);
		p.addBooleanArg("-cs", "Use score cache so omega files are only read once", false, true);
		p.addStringArg("-chr", "Single chromosome to use", false, null);
		p.parse(args);
		boolean debug = p.getBooleanArg("-d");
		if(debug) {
			logger.setLevel(Level.DEBUG);
			OmegaKmerFileReader.logger.setLevel(Level.DEBUG);
		}
		String omegaDirectory = p.getStringArg("-w");
		int kmerSize = p.getIntArg("-k");
		int fileChunkSize = p.getIntArg("-c");
		String geneBedFile = p.getStringArg("-g");
		double quantile = p.getDoubleArg("-q");
		String outQuantile = p.getStringArg("-oq");
		String outSigKmerBed = p.getStringArg("-os");
		double pvalCutoff = p.getDoubleArg("-p");
		boolean cache = p.getBooleanArg("-cs");
		String singleChr = p.getStringArg("-chr");
		
		KmerOmegaAnalysis ka = new KmerOmegaAnalysis(omegaDirectory, kmerSize, fileChunkSize, geneBedFile, cache);
		
		if(outQuantile != null) {
			ka.writeOmegaQuantiles(quantile, outQuantile, singleChr);
		}
		
		if(outSigKmerBed != null) {
			ka.writeMergedSignificantKmers(pvalCutoff, outSigKmerBed, singleChr);
		}
		
		logger.info("");
		logger.info("All done.");
	}

}
