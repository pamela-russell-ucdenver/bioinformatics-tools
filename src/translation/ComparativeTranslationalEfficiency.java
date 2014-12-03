package translation;

import java.io.FileWriter;
import java.io.IOException;

import org.apache.log4j.Logger;

import net.sf.samtools.util.CloseableIterator;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.AnnotationCollection;
import guttmanlab.core.util.CommandLineParser;
import guttmanlab.core.util.CountLogger;

public class ComparativeTranslationalEfficiency {
	
	private TranslationalEfficiency sample1;
	private TranslationalEfficiency sample2;
	private static Logger logger = Logger.getLogger(ComparativeTranslationalEfficiency.class.getName());
	
	public ComparativeTranslationalEfficiency(TranslationalEfficiency teSample1, TranslationalEfficiency teSample2) {
		sample1 = teSample1;
		sample2 = teSample2;
	}
	
	/**
	 * @param geneBed Bed file of genome annotation
	 * @param chrSizes Chromsome size file
	 * @param isStrandSpecific Whether the libraries are strand specific
	 * @param ribosomeBam1 Bam file of ribosome profiling sample 1
	 * @param controlBam1 Bam file of control sample 1
	 * @param ribosomeBam2 Bam file of ribosome profiling sample 2
	 * @param controlBam2 Bam file of control sample 2
	 * @param ribosomeGenomeTotal1 Optional total number of ribosome reads mapped to genome (instead of computing from data), sample 1
	 * @param controlGenomeTotal1 Optional total number of control reads mapped to genome (instead of computing from data), sample 1
	 * @param ribosomeGenomeTotal2 Optional total number of ribosome reads mapped to genome (instead of computing from data), sample 2
	 * @param controlGenomeTotal2 Optional total number of control reads mapped to genome (instead of computing from data), sample 2
	 * @param ribosomeExonTotal1 Optional total number of ribosome reads mapped to exons (instead of computing from data), sample 1
	 * @param controlExonTotal1 Optional total number of control reads mapped to exons (instead of computing from data), sample 1
	 * @param ribosomeExonTotal2 Optional total number of ribosome reads mapped to exons (instead of computing from data), sample 2
	 * @param controlExonTotal2 Optional total number of control reads mapped to exons (instead of computing from data), sample 2
	 * @throws IOException
	 */
	public static ComparativeTranslationalEfficiency factory(String geneBed, String chrSizes, boolean isStrandSpecific, 
			String ribosomeBam1, String controlBam1, String ribosomeBam2, String controlBam2, 
			double ribosomeGenomeTotal1, double controlGenomeTotal1, double ribosomeGenomeTotal2, double controlGenomeTotal2, 
			double ribosomeExonTotal1, double controlExonTotal1, double ribosomeExonTotal2, double controlExonTotal2) throws IOException {
		logger.info("");
		logger.info("Initializing sample 1...");
		TranslationalEfficiency te1 = new TranslationalEfficiency(ribosomeBam1, controlBam1, geneBed, chrSizes, ribosomeGenomeTotal1, controlGenomeTotal1, ribosomeExonTotal1, controlExonTotal1, isStrandSpecific);
		logger.info("");
		logger.info("Initializing sample 2...");
		TranslationalEfficiency te2 = new TranslationalEfficiency(ribosomeBam2, controlBam2, geneBed, chrSizes, ribosomeGenomeTotal2, controlGenomeTotal2, ribosomeExonTotal2, controlExonTotal2, isStrandSpecific);
		return new ComparativeTranslationalEfficiency(te1, te2);
	}
	
	/**
	 * Compute a ratio of TE scores provided that the gene is expressed in both control samples
	 * and at least one sample has a significant number of ribosome footprints in the CDS
	 * @param gene The gene
	 * @return Raw ratio of TE scores for the CDSs: sample 2 / sample 1
	 */
	public double getCdsTERawFoldChange(Gene gene) {
		// Check that at least one sample has ribosome footprints
		if(!sample1.cdsHasRibosomeFootprints(gene) && !sample2.cdsHasRibosomeFootprints(gene)) {
			return Double.NaN;
		}
		double te1 = sample1.getCdsTE(gene);
		if(te1 == Double.NaN) return Double.NaN;
		double te2 = sample2.getCdsTE(gene);
		if(te2 == Double.NaN) return Double.NaN;
		return te2 / te1;
	}
	
	/**
	 * @param gene The gene
	 * @param base Logarithm base
	 * @return Log ratio of TE scores for the CDSs: sample 2 / sample 1
	 */
	public double getCdsTELogFoldChange(Gene gene, double base) {
		double raw = getCdsTERawFoldChange(gene);
		return Math.log(raw) / Math.log(base);
	}
	
	/**
	 * Write TE scores for the CDS of each gene to a table
	 * @param geneBed Bed file of genes to analyze
	 * @param outputTable Table file to write
	 * @throws IOException 
	 */
	public void writeLogFoldChangeCdsTEsToTable(String geneBed, String chrSizeFile, String outputTable, double logBase) throws IOException {
		logger.info("");
		FileWriter w = new FileWriter(outputTable);
		AnnotationCollection<Gene> genes = BEDFileIO.loadFromFile(geneBed, chrSizeFile);
		int numGenes = genes.getNumAnnotations();
		CountLogger countLogger = new CountLogger(numGenes, 20);
		logger.info("Writing table with log fold change of CDS TE scores for " + numGenes + " genes in " + geneBed + "...");
		CloseableIterator<Gene> iter = genes.sortedIterator();
		String header = "gene\t";
		header += "TE_score_CDS_" + sample1.getSampleName() + "\t";
		header += "TE_score_CDS_" + sample2.getSampleName() + "\t";
		header += "log" + logBase + "_fold_change";
		w.write(header + "\n");
		while(iter.hasNext()) {
			Gene gene = iter.next();
			countLogger.advance();
			double te1 = sample1.getCdsTE(gene);
			double te2 = sample2.getCdsTE(gene);
			double logFoldChange = getCdsTELogFoldChange(gene, logBase);
			String line = gene.getName() + "\t";
			line += te1 + "\t";
			line += te2 + "\t";
			line += logFoldChange + "\t";
			w.write(line + "\n");
		}
		w.close();
		iter.close();
		logger.info("Done writing to " + outputTable + ".");
	}

	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-g", "Gene annotation bed file for computing totals", true);
		p.addStringArg("-c", "Chromosome size file", true);
		p.addStringArg("-gc", "Bed file of genes for computing log fold change TE of CDS", true);
		p.addStringArg("-ot", "Output table of log fold change CDS TE", true);
		p.addBooleanArg("-ss", "Libraries are strand specific", false, true);
		p.addDoubleArg("-lb", "Logarithm base", false, 2);
		p.addStringArg("-r1", "Ribosome bam file, sample 1", true);
		p.addStringArg("-m1", "Control bam file, sample 1", true);
		p.addDoubleArg("-mtg1", "Control global genome total (instead of computing from data), sample 1", false, -1);
		p.addDoubleArg("-rtg1", "Ribosome global genome total (instead of computing from data), sample 1", false, -1);
		p.addDoubleArg("-mte1", "Control global exon total (instead of computing from data), sample 1", false, -1);
		p.addDoubleArg("-rte1", "Ribosome global exon total (instead of computing from data), sample 1", false, -1);
		p.addStringArg("-r2", "Ribosome bam file, sample 2", true);
		p.addStringArg("-m2", "Control bam file, sample 2", true);
		p.addDoubleArg("-mtg2", "Control global genome total (instead of computing from data), sample 2", false, -1);
		p.addDoubleArg("-rtg2", "Ribosome global genome total (instead of computing from data), sample 2", false, -1);
		p.addDoubleArg("-mte2", "Control global exon total (instead of computing from data), sample 2", false, -1);
		p.addDoubleArg("-rte2", "Ribosome global exon total (instead of computing from data), sample 2", false, -1);
		p.parse(args);
		String geneAnnotationBed = p.getStringArg("-g");
		String chrSizes = p.getStringArg("-c");
		String geneBed = p.getStringArg("-gc");
		String outputTable = p.getStringArg("-ot");
		String ribosomeBam1 = p.getStringArg("-r1");
		String controlBam1 = p.getStringArg("-m1");
		double ribosomeGenomeTotal1 = p.getDoubleArg("-rtg1");
		double controlGenomeTotal1 = p.getDoubleArg("-mtg1");
		double ribosomeExonTotal1 = p.getDoubleArg("-rte1");
		double controlExonTotal1 = p.getDoubleArg("-mte1");
		boolean strandSpecific = p.getBooleanArg("-ss");
		String ribosomeBam2 = p.getStringArg("-r2");
		String controlBam2 = p.getStringArg("-m2");
		double ribosomeGenomeTotal2 = p.getDoubleArg("-rtg2");
		double controlGenomeTotal2 = p.getDoubleArg("-mtg2");
		double ribosomeExonTotal2 = p.getDoubleArg("-rte2");
		double controlExonTotal2 = p.getDoubleArg("-mte2");
		double logBase = p.getDoubleArg("-lb");
		
		ComparativeTranslationalEfficiency c = factory(geneAnnotationBed, chrSizes, strandSpecific, 
			ribosomeBam1, controlBam1, ribosomeBam2, controlBam2, 
			ribosomeGenomeTotal1, controlGenomeTotal1, ribosomeGenomeTotal2, controlGenomeTotal2, 
			ribosomeExonTotal1, controlExonTotal1, ribosomeExonTotal2, controlExonTotal2);
		
		c.writeLogFoldChangeCdsTEsToTable(geneBed, chrSizes, outputTable, logBase);
		
		logger.info("");
		logger.info("All done.");
		
	}

	
}
