package translation;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.apache.log4j.Logger;

import broad.core.math.ScanStatistics;
import net.sf.samtools.util.CloseableIterator;
import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.AnnotationCollection;
import guttmanlab.core.annotationcollection.BAMSingleReadCollection;
import guttmanlab.core.coordinatespace.CoordinateSpace;
import guttmanlab.core.util.CommandLineParser;
import guttmanlab.core.util.CountLogger;
import guttmanlab.core.util.StringParser;

public class TranslationalEfficiency {
	
	private double ribosomeGlobalGenomeTotal;
	private double controlGlobalGenomeTotal;
	private double ribosomeGlobalExonTotal;
	private double controlGlobalExonTotal;
	private double normalizationFactor; // Ribosome total / control total
	private BAMSingleReadCollection ribosomeData;
	private BAMSingleReadCollection controlData;
	private String chrSizeFile;
	private String controlBamFile;
	private String ribosomeBamFile;
	private long totalChrSize;
	private static final double EXPRESSION_SCAN_PVAL_CUTOFF = 0.05;
	private Map<String, Double> expressionScanPvals;
	private Map<String, Double> ribosomeCdsScanPvals;
	private Map<String, Double> ribosomeCounts;
	private Map<String, Double> controlCounts;
	private boolean strandSpecific;
	private double controlGlobalGenomeLambda;
	private double ribosomeGlobalGenomeLambda;
	private String sampleName;
	
	private static Logger logger = Logger.getLogger(TranslationalEfficiency.class.getName());
	
	/**
	 * @param ribosomeBam Bam file of ribosome profiling sample
	 * @param controlBam Bam file of control sample
	 * @param geneBed Bed file of genome annotation
	 * @param chrSizes Chromsome size file
	 * @param ribosomeGenomeTotal Optional total number of ribosome reads mapped to genome (instead of computing from data)
	 * @param controlGenomeTotal Optional total number of control reads mapped to genome (instead of computing from data)
	 * @param ribosomeExonTotal Optional total number of ribosome reads mapped to exons (instead of computing from data)
	 * @param controlExonTotal Optional total number of control reads mapped to exons (instead of computing from data)
	 * @param isStrandSpecific Whether the libraries are strand specific
	 * @throws IOException
	 */
	public TranslationalEfficiency(String ribosomeBam, String controlBam, String geneBed, String chrSizes, double ribosomeGenomeTotal, double controlGenomeTotal, double ribosomeExonTotal, double controlExonTotal, boolean isStrandSpecific) throws IOException {
		logger.info("");
		
		// Make sample name from ribosome bam file
		StringParser s = new StringParser();
		s.parse(ribosomeBam, "/");
		String withExtension = s.asString(s.getFieldCount() - 1);
		String s1 = withExtension.replaceAll(".sorted.bam", "");
		sampleName = s1.replaceAll(".bam", "");
		logger.info("Sample name is " + sampleName + ".");
		
		// Set whether the libraries are strand specific (affects counts over annotations)
		strandSpecific = isStrandSpecific;
		
		// Initialize caches of expression P values and region counts
		expressionScanPvals = new HashMap<String, Double>();
		ribosomeCdsScanPvals = new HashMap<String, Double>();
		ribosomeCounts = new HashMap<String, Double>();
		controlCounts = new HashMap<String, Double>();
		
		// Save chromosome size file to use when loading annotations
		chrSizeFile = chrSizes;
		
		// Load read mapping data
		controlBamFile = controlBam;
		ribosomeBamFile = ribosomeBam;
		ribosomeData = new BAMSingleReadCollection(new File(ribosomeBamFile));
		controlData = new BAMSingleReadCollection(new File(controlBamFile));

		// Compute global read counts
		controlGlobalGenomeTotal = controlGenomeTotal;
		if(controlGlobalGenomeTotal <= 0) {
			logger.info("Computing total genome read count for control sample...");
			controlGlobalGenomeTotal = controlData.getNumAnnotations();
		}
		logger.info(controlGlobalGenomeTotal + " total reads in control fraction.");
		ribosomeGlobalGenomeTotal = ribosomeGenomeTotal;
		if(ribosomeGlobalGenomeTotal <= 0) {
			logger.info("Computing total genome read count for ribosome sample...");
			ribosomeGlobalGenomeTotal = ribosomeData.getNumAnnotations();
		}
		logger.info(ribosomeGlobalGenomeTotal + " total reads in ribosome fraction.");

		// Compute scan distribution parameters
		CoordinateSpace chrs = new CoordinateSpace(chrSizeFile);
		totalChrSize = chrs.getTotalReferenceLength();
		controlGlobalGenomeLambda = controlGlobalGenomeTotal / totalChrSize;
		logger.info("Total chromosome size is " + totalChrSize + ". Global mapped reads in control fraction is " + controlGlobalGenomeTotal + ". Control global lambda is " + controlGlobalGenomeLambda + ".");
		ribosomeGlobalGenomeLambda = ribosomeGlobalGenomeTotal / totalChrSize;
		
		// Calculate TE normalization factor
		AnnotationCollection<Gene> genes = BEDFileIO.loadFromFile(geneBed, chrSizeFile);
		if(controlExonTotal < 0) {
			logger.info("Computing total exon read count for control sample...");
			controlExonTotal = exonTotalcontrol(genes);
		}
		controlGlobalExonTotal = controlExonTotal;
		logger.info(controlGlobalExonTotal + " total exon reads in control fraction.");
		if(ribosomeExonTotal < 0) {
			logger.info("Computing total exon read count for ribosome sample...");
			ribosomeExonTotal = exonTotalRibosome(genes);
		}
		ribosomeGlobalExonTotal = ribosomeExonTotal;
		logger.info(ribosomeGlobalExonTotal + " total exon reads in ribosome fraction.");
		normalizationFactor = ribosomeGlobalExonTotal / controlGlobalExonTotal;
	}
	
	public String getSampleName() {
		return sampleName;
	}
	
	/**
	 * @param ribosomeBam Bam file of ribosome profiling sample
	 * @param controlBam Bam file of control sample
	 * @param geneBed Bed file of genome annotation
	 * @param chrSizes Chromsome size file
	 * @param ribosomeGenomeTotal Optional total number of ribosome reads mapped to genome (instead of computing from data)
	 * @param controlGenomeTotal Optional total number of control reads mapped to genome (instead of computing from data)
	 * @param ribosomeExonTotal Optional total number of ribosome reads mapped to exons (instead of computing from data)
	 * @param controlExonTotal Optional total number of control reads mapped to exons (instead of computing from data)
	 * @param isStrandSpecific Whether the libraries are strand specific
	 * @throws IOException
	 */
	public static TranslationalEfficiency factory(String ribosomeBam, String controlBam, String geneBed, String chrSizes, double ribosomeGenomeTotal, double controlGenomeTotal, double ribosomeExonTotal, double controlExonTotal, boolean isStrandSpecific) throws IOException {
		return new TranslationalEfficiency(ribosomeBam, controlBam, geneBed, chrSizes, ribosomeGenomeTotal, controlGenomeTotal, ribosomeExonTotal, controlExonTotal, isStrandSpecific);
	}
	
	/**
	 * @param gene Gene
	 * @return Whether the gene contains a significant number of control reads compared to genomic background
	 */
	public boolean isExpressed(Annotation gene) {
		double expressionPval = getExpressionScanPval(gene);
		return expressionPval < EXPRESSION_SCAN_PVAL_CUTOFF;
	}
	
	/**
	 * @param gene Gene
	 * @return Whether the CDS of the gene contains a significant number of ribosome footprints compared to genomic background
	 */
	public boolean cdsHasRibosomeFootprints(Gene gene) {
		double ribosomePval = getCdsRibosomeScanPval(gene);
		return ribosomePval < EXPRESSION_SCAN_PVAL_CUTOFF;
	}
	
	/**
	 * Get translational efficiency score for a region
	 * @param region The region
	 * @return The TE score
	 */
	public double getTE(Annotation region, Annotation parentGene) {
		if(!isExpressed(parentGene)) {
			return Double.NaN;
		}
		double ribosomeReads = getRibosomeCount(region);
		double controlReads = getControlCount(region);
		return (ribosomeReads / controlReads) / normalizationFactor;
	}
	
	/**
	 * Get translational efficiency score for the CDS of a gene
	 * @param gene The gene
	 * @return The TE score of the CDS
	 */
	public double getCdsTE(Gene gene) {
		Annotation cds = gene.getCodingRegion();
		return getTE(cds, gene);
	}
	
	/**
	 * Write TE scores for the CDS of each gene to a bed file
	 * @param geneBed Bed file of genes to analyze
	 * @param outputBed Bed file to write
	 * @throws IOException 
	 */
	public void writeCdsTEsToBed(String geneBed, String outputBed) throws IOException {
		logger.info("");
		FileWriter w = new FileWriter(outputBed);
		AnnotationCollection<Gene> genes = BEDFileIO.loadFromFile(geneBed, chrSizeFile);
		int numGenes = genes.getNumAnnotations();
		CountLogger countLogger = new CountLogger(numGenes, 20);
		logger.info("Writing bed file with CDS TE scores for " + genes.getNumAnnotations() + " genes in " + geneBed + "...");
		CloseableIterator<Gene> iter = genes.sortedIterator();
		while(iter.hasNext()) {
			countLogger.advance();
			Gene gene = iter.next();
			double te = getCdsTE(gene);
			String bedLine = gene.toBED(te);
			w.write(bedLine + "\n");
		}
		w.close();
		iter.close();
		logger.info("Done writing to " + outputBed + ".");
	}
	
	
	/**
	 * Write TE scores for the CDS of each gene to a table
	 * @param geneBed Bed file of genes to analyze
	 * @param outputTable Table file to write
	 * @throws IOException 
	 */
	public void writeCdsTEsToTable(String geneBed, String outputTable) throws IOException {
		logger.info("");
		FileWriter w = new FileWriter(outputTable);
		AnnotationCollection<Gene> genes = BEDFileIO.loadFromFile(geneBed, chrSizeFile);
		int numGenes = genes.getNumAnnotations();
		CountLogger countLogger = new CountLogger(numGenes, 20);
		logger.info("Writing table with CDS TE scores for " + numGenes + " genes in " + geneBed + "...");
		CloseableIterator<Gene> iter = genes.sortedIterator();
		String header = "gene\t";
		header += "global_genome_count_control\t";
		header += "global_genome_count_ribosome\t";
		header += "global_exon_count_control\t";
		header += "global_exon_count_ribosome\t";
		header += "global_lambda\t";
		header += "gene_size\t";
		header += "gene_lambda\t";
		header += "expression_scan_pval\t";
		header += "gene_count_control\t";
		header += "gene_count_ribosome\t";
		header += "cds_count_control\t";
		header += "cds_count_ribosome\t";
		header += "TE_normalization_factor\t";
		header += "TE_score_CDS\t";
		w.write(header + "\n");
		while(iter.hasNext()) {
			Gene gene = iter.next();
			countLogger.advance();
			Annotation cds = gene.getCodingRegion();
			double expPval = getExpressionScanPval(gene);
			double rCountCDS = getRibosomeCount(cds);
			double mCountCDS = getControlCount(cds);
			double mCountGene = getControlCount(gene);
			double rCountGene = getRibosomeCount(gene);
			double te = getCdsTE(gene);
			int geneSize = gene.size();
			double geneLambda = mCountGene / geneSize;
			String line = gene.getName() + "\t";
			line += controlGlobalGenomeTotal + "\t";
			line += ribosomeGlobalGenomeTotal + "\t";
			line += controlGlobalExonTotal + "\t";
			line += ribosomeGlobalExonTotal + "\t";
			line += controlGlobalGenomeLambda + "\t";
			line += geneSize + "\t";
			line += geneLambda + "\t";
			line += expPval + "\t";
			line += mCountGene + "\t";
			line += rCountGene + "\t";
			line += mCountCDS + "\t";
			line += rCountCDS + "\t";
			line += normalizationFactor + "\t";
			line += te + "\t";
			w.write(line + "\n");
		}
		w.close();
		iter.close();
		logger.info("Done writing to " + outputTable + ".");
	}
	
	/**
	 * @param reads Read mappings
	 * @param genes Gene annotation
	 * @return Total number of reads mapped to exons
	 */
	private double exonTotalRibosome(AnnotationCollection<Gene> genes) {
		int numGenes = genes.getNumAnnotations();
		CountLogger countLogger = new CountLogger(numGenes, 20);
		CloseableIterator<Gene> iter = genes.sortedIterator();
		double total = 0;
		while(iter.hasNext()) {
			Gene gene = iter.next();
			total += getRibosomeCount(gene);
			countLogger.advance();
		} 
		iter.close();
		return total;
	}
	
	/**
	 * @param reads Read mappings
	 * @param genes Gene annotation
	 * @return Total number of reads mapped to exons
	 */
	private double exonTotalcontrol(AnnotationCollection<Gene> genes) {
		int numGenes = genes.getNumAnnotations();
		CountLogger countLogger = new CountLogger(numGenes, 20);
		CloseableIterator<Gene> iter = genes.sortedIterator();
		double total = 0;
		while(iter.hasNext()) {
			Gene gene = iter.next();
			total += getControlCount(gene);
			countLogger.advance();
		} 
		iter.close();
		return total;
	}
	
	/**
	 * @param gene Gene
	 * @return Scan P value for read count over gene in control sample
	 */
	public double getExpressionScanPval(Annotation gene) {
		if(expressionScanPvals.containsKey(gene.toBED())) {
			return expressionScanPvals.get(gene.toBED()).doubleValue();
		}
		int controlCount = (int)getControlCount(gene);
		int geneSize = gene.size();
		double rtrn = ScanStatistics.calculatePVal(controlCount, controlGlobalGenomeLambda, geneSize, totalChrSize);
		expressionScanPvals.put(gene.toBED(), Double.valueOf(rtrn));
		return rtrn;
	}
	
	
	/**
	 * @param gene Gene
	 * @return Scan P value for ribosome footprint count over CDS of gene
	 */
	public double getCdsRibosomeScanPval(Gene gene) {
		if(ribosomeCdsScanPvals.containsKey(gene.toBED())) {
			return ribosomeCdsScanPvals.get(gene.toBED()).doubleValue();
		}
		Annotation cds = gene.getCodingRegion();
		int ribosomeCdsCount = (int)getRibosomeCount(cds);
		int cdsSize = cds.size();
		double rtrn = ScanStatistics.calculatePVal(ribosomeCdsCount, ribosomeGlobalGenomeLambda, cdsSize, totalChrSize);
		ribosomeCdsScanPvals.put(gene.toBED(), Double.valueOf(rtrn));
		return rtrn;
	}
	
	
	/**
	 * @param gene Gene
	 * @return Ribosome read count over gene
	 */
	public double getRibosomeCount(Annotation gene) {
		Gene geneToUse = new Gene(gene);
		if(!strandSpecific) {
			geneToUse.setOrientation(Strand.BOTH);
		}
		if(ribosomeCounts.containsKey(geneToUse.toBED())) {
			return ribosomeCounts.get(geneToUse.toBED()).doubleValue();
		}
		int ribosome = ribosomeData.numOverlappers(geneToUse, false);
		ribosomeCounts.put(geneToUse.toBED(), Double.valueOf(ribosome));
		return ribosome;
	}
	
	/**
	 * @param gene Gene
	 * @return control read count over gene
	 */
	public double getControlCount(Annotation gene) {
		Gene geneToUse = new Gene(gene);
		if(!strandSpecific) {
			geneToUse.setOrientation(Strand.BOTH);
		}
		if(controlCounts.containsKey(geneToUse.toBED())) {
			return controlCounts.get(geneToUse.toBED()).doubleValue();
		}
		int rtrn = controlData.numOverlappers(geneToUse, false);
		controlCounts.put(geneToUse.toBED(), Double.valueOf(rtrn));
		return rtrn;
	}
	
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-r", "Ribosome bam file", true);
		p.addStringArg("-m", "Control bam file", true);
		p.addStringArg("-g", "Gene annotation bed file for computing totals", true);
		p.addStringArg("-c", "Chromosome size file", true);
		p.addStringArg("-gc", "Bed file of genes for computing TE of CDS", true);
		p.addStringArg("-ob", "Output bed file with gene score set to CDS TE", false, null);
		p.addStringArg("-ot", "Output table of CDS TE", false, null);
		p.addDoubleArg("-mtg", "Control global genome total (instead of computing from data)", false, -1);
		p.addDoubleArg("-rtg", "Ribosome global genome total (instead of computing from data)", false, -1);
		p.addDoubleArg("-mte", "Control global exon total (instead of computing from data)", false, -1);
		p.addDoubleArg("-rte", "Ribosome global exon total (instead of computing from data)", false, -1);
		p.addBooleanArg("-ss", "Libraries are strand specific", false, true);
		p.parse(args);
		String ribosomeBam = p.getStringArg("-r");
		String controlBam = p.getStringArg("-m");
		String geneAnnotationBed = p.getStringArg("-g");
		String chrSizes = p.getStringArg("-c");
		String geneBed = p.getStringArg("-gc");
		String outputBed = p.getStringArg("-ob");
		String outputTable = p.getStringArg("-ot");
		double ribosomeGenomeTotal = p.getDoubleArg("-rtg");
		double controlGenomeTotal = p.getDoubleArg("-mtg");
		double ribosomeExonTotal = p.getDoubleArg("-rte");
		double controlExonTotal = p.getDoubleArg("-mte");
		boolean strandSpecific = p.getBooleanArg("-ss");
		
		TranslationalEfficiency te = new TranslationalEfficiency(ribosomeBam, controlBam, geneAnnotationBed, chrSizes, ribosomeGenomeTotal, controlGenomeTotal, ribosomeExonTotal, controlExonTotal, strandSpecific);
		
		if(outputBed != null) te.writeCdsTEsToBed(geneBed, outputBed);
		
		if(outputTable != null) te.writeCdsTEsToTable(geneBed, outputTable);
		
		logger.info("");
		logger.info("All done.");
		
	}
	
}
