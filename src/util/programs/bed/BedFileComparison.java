package util.programs.bed;

import guttmanlab.core.util.CommandLineParser;

import java.io.File;
import java.io.IOException;
import java.util.Collection;

import nextgen.core.annotation.Gene;

import org.apache.log4j.Logger;

import broad.pda.annotation.BEDFileParser;

public class BedFileComparison {
	
	private static Logger logger = Logger.getLogger(BedFileComparison.class.getName());
	
	private static void eraseInfoExceptGeneModel(Gene gene) {
		gene.setBedScore(0);
		gene.setCountScore(0);
		gene.setName("");
		gene.setScore(0);
	}
	
	public static void compareGeneModels(String bed1, String bed2) throws IOException {
		Collection<Gene> genes1 = BEDFileParser.loadData(new File(bed1));
		Collection<Gene> genes2 = BEDFileParser.loadData(new File(bed2));
		for(Gene gene : genes1) {
			eraseInfoExceptGeneModel(gene);
		}
		for(Gene gene : genes2) {
			eraseInfoExceptGeneModel(gene);
		}
		int bed1only = 0;
		int bed2only = 0;
		int shared = 0;
		for(Gene gene : genes1) {
			if(genes2.contains(gene)) {
				shared++;
				logger.info("Shared\t" + gene.toUCSC());
			} else {
				bed1only++;
			}
		}
		for(Gene gene : genes2) {
			if(!genes1.contains(gene)) {
				bed2only++;
			}
		}
		logger.info(bed1 + " only\t" + bed1only);
		logger.info(bed2 + " only\t" + bed2only);
		logger.info("Shared\t" + shared);
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-b1", "Bed file 1", true);
		p.addStringArg("-b2", "Bed file 2", true);
		p.parse(args);
		String bed1 = p.getStringArg("-b1");
		String bed2 = p.getStringArg("-b2");
		
		compareGeneModels(bed1, bed2);
		
	}

}
