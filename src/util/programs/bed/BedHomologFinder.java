package util.programs.bed;

import guttmanlab.core.util.CommandLineParser;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeSet;

import org.apache.log4j.Logger;

import broad.pda.annotation.BEDFileParser;

import nextgen.core.annotation.Gene;

public class BedHomologFinder {
	
	private Map<String, Collection<Gene>> genes;
	private static Logger logger = Logger.getLogger(BedHomologFinder.class.getName());
	
	private BedHomologFinder(String annotationBedFile) throws IOException {
		genes = BEDFileParser.loadDataByChr(new File(annotationBedFile));
	}
	
	private Collection<Gene> findHomologs(Gene gene, double pctExonOverlap) {
		if(pctExonOverlap <= 0 || pctExonOverlap > 1) {
			throw new IllegalArgumentException("Minimum percent overlap must be between 0 and 1");
		}
		Collection<Gene> rtrn = new TreeSet<Gene>();
		for(Gene other : genes.get(gene.getChr())) {
			if(gene.percentOverlapping(other) >= pctExonOverlap && gene.getOrientation().equals(other.getOrientation())) {
				rtrn.add(other);
			}
		}
		return rtrn;
	}
	
	private void writeHomologs(Gene gene, double pctExonOverlap, FileWriter w) throws IOException {
		Collection<Gene> homologs = findHomologs(gene, pctExonOverlap);
		for(Gene homolog : homologs) {
			w.write(gene.getName() + "\t" + homolog.getName() + "\n");
		}
	}
	
	private void writeHomologs(String genesBed, double pctExonOverlap, String outFile) throws IOException {
		logger.info("");
		logger.info("Writing table of homologs to " + outFile);
		Map<String, Collection<Gene>> geneSet = BEDFileParser.loadDataByChr(new File(genesBed));
		FileWriter w = new FileWriter(outFile);
		for(String chr : geneSet.keySet()) {
			logger.info(chr);
			for(Gene gene : geneSet.get(chr)) {
				writeHomologs(gene, pctExonOverlap, w);
			}
		}
		w.close();
		logger.info("Done writing table.");
	}
	
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-ba", "Bed gene annotation", true);
		p.addStringArg("-bg", "Bed file of genes to find homologs for", true);
		p.addStringArg("-o", "Output table", true);
		p.addDoubleArg("-mp", "Minimum percent exon overlap", true);
		p.parse(args);
		String bedAnnotation = p.getStringArg("-ba");
		String bedGenes = p.getStringArg("-bg");
		String outTable = p.getStringArg("-o");
		double pctOverlap = p.getDoubleArg("-mp");
		
		BedHomologFinder b = new BedHomologFinder(bedAnnotation);
		b.writeHomologs(bedGenes, pctOverlap, outTable);
		
		logger.info("");
		logger.info("All done.");
		
	}
	
}
