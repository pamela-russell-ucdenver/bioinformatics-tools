package util.programs.bed;

import guttmanlab.core.util.CommandLineParser;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;

import broad.pda.annotation.BEDFileParser;
import nextgen.core.annotation.Gene;
import nextgen.core.annotation.Annotation;

/**
 * A special class to write the target region for CRISPR V6.5 experiment July 2014
 * @author prussell
 *
 */
public class CrisprCDSWriter {

	/**
	 * @param gene Gene
	 * @return The CDS or the full gene if untranslated
	 */
	private static Gene getCDSIfTranslated(Gene gene) {
		Gene cds = gene.getCDS();
		if(cds != null) return cds;
		return gene;
	}
	
	/**
	 * @param gene Gene
	 * @return New gene with last exon removed or original gene if only one exon
	 */
	private static Gene removeLastExonIfMultiExon(Gene gene) {
		if(gene.numBlocks() < 2) return gene;
		Collection<? extends Annotation> exons = gene.getExonSet();
		exons.remove(gene.get3PrimeExon());
		return new Gene(exons, gene.getName());
	}

	/**
	 * @param gene
	 * @return The transformation of the gene requested by Patrick, July 2014
	 */
	private static Gene transformGene(Gene gene) {
		Gene cds = getCDSIfTranslated(gene);
		Gene noLastExon = removeLastExonIfMultiExon(cds);
		Gene rtrn = noLastExon;
		rtrn.setName(gene.getName());
		return rtrn;
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-b", "Input bed", true);
		p.addStringArg("-o", "Output bed", true);
		p.parse(args);
		String input = p.getStringArg("-b");
		String output = p.getStringArg("-o");
		
		Map<String, Collection<Gene>> genes = BEDFileParser.loadDataByChr(input);
		FileWriter w = new FileWriter(output);
		for(String chr : genes.keySet()) {
			for(Gene gene : genes.get(chr)) {
				w.write(transformGene(gene).toBED() + "\n");
			}
		}
		w.close();
		
	}

}
