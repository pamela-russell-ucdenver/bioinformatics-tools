/**
 * 
 */
package util.programs.counts;

import guttmanlab.core.util.CommandLineParser;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;



import org.apache.log4j.Logger;

import broad.core.math.Statistics;
import broad.pda.annotation.BEDFileParser;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Gene;
import nextgen.core.coordinatesystem.GenomicSpace;
import nextgen.core.coordinatesystem.TranscriptomeSpace;
import nextgen.core.feature.GeneWindow;
import nextgen.core.model.ScanStatisticDataAlignmentModel;
import nextgen.core.model.score.ScanStatisticScore;
import nextgen.core.readFilters.GenomicSpanFilter;
import nextgen.core.readFilters.ProperPairFilter;

/**
 * @author prussell
 * Compare an RNA seq sample to background expression
 * Compare read count ratio in individual features to overall ratio
 */
public class EnrichmentByFeature {
	
	private ScanStatisticDataAlignmentModel backgroundTranscriptomeData;
	private ScanStatisticDataAlignmentModel signalTranscriptomeData;
	private ScanStatisticDataAlignmentModel backgroundGenomicData;
	private ScanStatisticDataAlignmentModel signalGenomicData;
	private Collection<Gene> genes;
	static Logger logger = Logger.getLogger(EnrichmentByFeature.class.getName());
	private TranscriptomeSpace transcriptomeSpace;
	private GenomicSpace genomeSpace;
	private static int DEFAULT_MAX_GENOMIC_SPAN = 300000;
	private static boolean DEFAULT_FULLY_CONTAINED = false;
	private boolean fullyContained;
	
	private EnrichmentByFeature(String backgroundBamFile, String signalBamFile, String bedFile, String chrSizeFile) throws IOException {
		this(backgroundBamFile, signalBamFile, bedFile, chrSizeFile, DEFAULT_MAX_GENOMIC_SPAN, DEFAULT_FULLY_CONTAINED);
	}
	
	private EnrichmentByFeature(String backgroundBamFile, String signalBamFile, String bedFile, String chrSizeFile, int maxGenomicSpan, boolean fully_contained) throws IOException {
		fullyContained = fully_contained;
		genomeSpace = new GenomicSpace(chrSizeFile);
		transcriptomeSpace = new TranscriptomeSpace(BEDFileParser.loadDataByChr(new File (bedFile)));
		backgroundTranscriptomeData = new ScanStatisticDataAlignmentModel(backgroundBamFile, transcriptomeSpace);
		backgroundGenomicData = new ScanStatisticDataAlignmentModel(backgroundBamFile, genomeSpace);
		signalGenomicData = new ScanStatisticDataAlignmentModel(signalBamFile, genomeSpace);
		signalTranscriptomeData = new ScanStatisticDataAlignmentModel(signalBamFile, transcriptomeSpace);		
		genes = BEDFileParser.loadData(new File(bedFile));
		// Add read filters
		backgroundTranscriptomeData.addFilter(new ProperPairFilter());
		backgroundTranscriptomeData.addFilter(new GenomicSpanFilter(maxGenomicSpan));
		signalTranscriptomeData.addFilter(new ProperPairFilter());
		signalTranscriptomeData.addFilter(new GenomicSpanFilter(maxGenomicSpan));
		
	}
	
	/**
	 * Calculate enrichment off of counts
	 * @param window 
	 * @param useGenomicSpace 
	 * @return The ratio divided by global count ratio
	 */
	public double calculateEnrichment(Annotation window, boolean useGenomicSpace) {
		
		double windowSignalAverageCoverage = 0;
		double windowBackgroundAverageCoverage = 0;
		
		if(!useGenomicSpace) {
			ScanStatisticScore s = new ScanStatisticScore(signalTranscriptomeData, window, signalTranscriptomeData.getCount(window), window.getSize(), false);
			windowSignalAverageCoverage = s.getAverageCoverage(signalTranscriptomeData);
			ScanStatisticScore b = new ScanStatisticScore(backgroundTranscriptomeData, window, backgroundTranscriptomeData.getCount(window), window.getSize(), false);
			windowBackgroundAverageCoverage = b.getAverageCoverage(backgroundTranscriptomeData);
		} else {
			ScanStatisticScore s = new ScanStatisticScore(signalGenomicData, window, signalGenomicData.getCount(window), window.getSize(), false);
			windowSignalAverageCoverage = s.getAverageCoverage(signalGenomicData);
			ScanStatisticScore b = new ScanStatisticScore(backgroundGenomicData, window, backgroundGenomicData.getCount(window), window.getSize(), false);
			windowBackgroundAverageCoverage = b.getAverageCoverage(backgroundGenomicData);
		}
		
		return (windowSignalAverageCoverage / windowBackgroundAverageCoverage);
	}

	/**
	 * Compute the max enrichment by windows in each feature
	 * Calculate the median for all genes for each feature type
	 * Write the median to a file
	 * @param outFile Output file
	 * @param windowSize Window size
	 * @param stepSize Step size
	 * @throws IOException
	 */
	private void computeAndWriteMaxWindowEnrichments(String outFile, int windowSize, int stepSize) throws IOException {
		FileWriter w = new FileWriter(outFile);
		logger.info("Calculating enrichments for " + genes.size() + " genes...");
		int numDone = 0;
		int numExpressed = 0;
		int numCoding = 0;
		int numNoncoding = 0;
		
		ArrayList<Double> codingTranscriptMaxWindowEnrichments = new ArrayList<Double>();
		ArrayList<Double> noncodingTranscriptMaxWindowEnrichments = new ArrayList<Double>();
		ArrayList<Double> cdsMaxWindowEnrichments = new ArrayList<Double>();
		ArrayList<Double> utr5MaxWindowEnrichments = new ArrayList<Double>();
		ArrayList<Double> utr3MaxWindowEnrichments = new ArrayList<Double>();
		ArrayList<Double> intronMaxWindowEnrichments = new ArrayList<Double>();
		
		for(Gene gene : genes) {

			numDone++;
			if(numDone % 1000 == 0) {
				logger.info("Finished " + numDone + " genes of which " + numExpressed + " are significantly expressed.");
			}
			
			if(!backgroundTranscriptomeData.isExpressed(gene, 0.05)) {
				continue;
			}
			
			numExpressed++;

			Gene cds = gene.getCDS();
			Gene utr5 = gene.get5UTRGene();
			Gene utr3 = gene.get3UTRGene();
			Gene introns = gene.getIntrons();					
			
			String info = "";
			
			boolean coding = gene.isCodingGene();
			if(coding) {
				info += "CODING\t" + gene.getName() + "\t";
				numCoding++;
			}
			else {
				info += "NONCODING\t" + gene.getName() + "\t";
				numNoncoding++;
			}
			
			// Entire transcript
			try {
				double maxEnrich = getMaxWindowEnrichment(gene, windowSize, stepSize);
				if(coding) codingTranscriptMaxWindowEnrichments.add(Double.valueOf(maxEnrich));
				else noncodingTranscriptMaxWindowEnrichments.add(Double.valueOf(maxEnrich));
				info += "transcript=" + maxEnrich + "\t";
			} catch(Exception e) {
				// 
			}
			
			
			if(coding) {
					// CDS
					try {
						double maxEnrich = getMaxWindowEnrichment(cds, windowSize, stepSize);
						cdsMaxWindowEnrichments.add(Double.valueOf(maxEnrich));
						info += "cds=" + maxEnrich + "\t";
					} catch(Exception e) {
						// 
					}
					// 5'UTR
					try {
						double maxEnrich = getMaxWindowEnrichment(utr5, windowSize, stepSize);
						utr5MaxWindowEnrichments.add(Double.valueOf(maxEnrich));
						info += "utr5=" + maxEnrich + "\t";
					} catch(Exception e) {
						// 
					}
					// 3'UTR
					try {
						double maxEnrich = getMaxWindowEnrichment(utr3, windowSize, stepSize);
						utr3MaxWindowEnrichments.add(Double.valueOf(maxEnrich));
						info += "utr3=" + maxEnrich + "\t";
					} catch(Exception e) {
						// 
					}
			}
			
			// Introns
			if(introns != null) {
				double maxEnrich = getMaxWindowEnrichment(introns, windowSize, stepSize, true);
				intronMaxWindowEnrichments.add(Double.valueOf(maxEnrich));
				info += "introns=" + maxEnrich + "\t";
			}
			
			logger.info(info);

		}
		
		w.write("Features\tMedianOfMaxWindowEnrichment\n");
		w.write("CodingTranscripts\t" + Statistics.median(codingTranscriptMaxWindowEnrichments) + "\n");
		w.write("NoncodingTranscripts\t" + Statistics.median(noncodingTranscriptMaxWindowEnrichments) + "\n");
		w.write("CDSs\t" + Statistics.median(cdsMaxWindowEnrichments) + "\n");
		w.write("5'UTRs\t" + Statistics.median(utr5MaxWindowEnrichments) + "\n");
		w.write("3'UTRs\t" + Statistics.median(utr3MaxWindowEnrichments) + "\n");
		w.write("Introns\t" + Statistics.median(intronMaxWindowEnrichments) + "\n");
		logger.info("Wrote median of max window enrichments to file " + outFile + ".");
		logger.info("There were " + numCoding + " expressed coding genes and " + numNoncoding + " expressed noncoding genes.");

		w.close();
	}
	
	/**
	 * Get the max enrichment over all windows in the gene
	 * @param gene
	 * @param windowSize
	 * @param stepSize
	 * @return
	 * @throws IOException
	 */
	private double getMaxWindowEnrichment(Gene gene, int windowSize, int stepSize) throws IOException {
		return getMaxWindowEnrichment(gene, windowSize, stepSize, false);
	}
	
	/**
	 * Get the max enrichment over all windows in the gene (count in either transcriptome space or genomic space)
	 * @param gene The gene
	 * @param windowSize Window size
	 * @param stepSize Step size
	 * @param useGenomicSpace Get counts in genomic space
	 * @return The enrichment of the max enriched window
	 * @throws IOException
	 */
	private double getMaxWindowEnrichment(Gene gene, int windowSize, int stepSize, boolean useGenomicSpace) throws IOException {
		Collection<GeneWindow> windows = gene.getWindows(windowSize, stepSize, 0);
		if(windows.size() == 0) throw new IllegalArgumentException("No windows");
		double rtrn = 0;
		for(GeneWindow window : windows) {
			double backgroundCount = 0;
			if(!useGenomicSpace) backgroundCount = backgroundTranscriptomeData.getCount(window, fullyContained);
			else {
				for(Annotation exon : window.getExonSet()) {
					backgroundCount += backgroundGenomicData.getCount(exon, fullyContained);
				}
			}
			if(backgroundCount == 0) continue;
			double enrichment = calculateEnrichment(window, useGenomicSpace);
			if(enrichment > rtrn) rtrn = enrichment;
		}
		return rtrn;
	}


	private void computeAndWriteEnrichments(String outMedianFile, String outTableFile) throws IOException {
		FileWriter medianWriter = new FileWriter(outMedianFile);
		FileWriter tableWriter = new FileWriter(outTableFile);
		
		int numExpressed = 0;
		
		logger.info("Calculating enrichments for " + genes.size() + " genes...");
		int numDone = 0;

		ArrayList<Double> cdsEnrichments = new ArrayList<Double>();
		ArrayList<Double> utr5Enrichments = new ArrayList<Double>();
		ArrayList<Double> utr3Enrichments = new ArrayList<Double>();
		ArrayList<Double> intronEnrichments = new ArrayList<Double>();
		
		for(Gene gene : genes) {

			numDone++;
			if(numDone % 1000 == 0) {
				logger.info("Finished " + numDone + " genes of which " + numExpressed + " are significantly expressed.");
			}
			
			if(!backgroundTranscriptomeData.isExpressed(gene, 0.05)) {
				continue;
			}
			
			numExpressed++;

			Gene cds = gene.getCDS();
			Gene utr5 = gene.get5UTRGene();
			Gene utr3 = gene.get3UTRGene();
			Gene introns = gene.getIntrons();
			
			
			
			if(cds != null) {
				if(cds.getSize() > 2) {
					// Gene is coding
					try {
						double enrichment = calculateEnrichment(cds, false);
						cdsEnrichments.add(Double.valueOf(enrichment));
						tableWriter.write("CDS\t" + cds.getName() + "\t" + enrichment + "\n");
					} catch(Exception e) {
						logger.warn("Caught exception. Skipping CDS for gene " + gene.getName());
					}
					try {
						double enrichment = calculateEnrichment(utr5, false);
						utr5Enrichments.add(Double.valueOf(enrichment));
						tableWriter.write("5UTR\t" + utr5.getName() + "\t" + enrichment + "\n");
					} catch(Exception e) {
						logger.warn("Caught exception. Skipping 5'UTR for gene " + gene.getName());
					}
					try {
						double enrichment = calculateEnrichment(utr3, false);
						utr3Enrichments.add(Double.valueOf(enrichment));
						tableWriter.write("3UTR\t" + utr3.getName() + "\t" + enrichment + "\n");
					} catch(Exception e) {
						logger.warn("Caught exception. Skipping 3'UTR for gene " + gene.getName());
					}
				}
			}
			
			if(introns != null) {
				try {
					double enrichment = calculateEnrichment(introns, false);
					intronEnrichments.add(Double.valueOf(enrichment));
					tableWriter.write("INTRONS\t" + introns.getName() + "\t" + enrichment + "\n");
				} catch(Exception e) {
					logger.warn("Caught exception. Skipping introns for gene " + gene.getName());
				}
			}
			

		}
			
		double medianCdsEnrichment = Statistics.median(cdsEnrichments);
		double medianUtr5Enrichment = Statistics.median(utr5Enrichments);
		double medianUtr3Enrichment = Statistics.median(utr3Enrichments);
		double medianIntronEnrichment = Statistics.median(intronEnrichments);
		
		medianWriter.write("Features\tMedianEnrichment\n");
		medianWriter.write("CDS\t" + medianCdsEnrichment + "\n");
		medianWriter.write("5UTR\t" + medianUtr5Enrichment + "\n");
		medianWriter.write("3UTR\t" + medianUtr3Enrichment + "\n");
		medianWriter.write("Introns\t" + medianIntronEnrichment + "\n");

		logger.info("Wrote median enrichments to file " + outMedianFile + ".");
		logger.info("Wrote table of feature enrichments to file " + outTableFile + ".");
		medianWriter.close();
		tableWriter.close();
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-b", "Background bam file", true);
		p.addStringArg("-s", "Signal bam file", true);
		p.addStringArg("-g", "Genes bed file", true);
		p.addStringArg("-c", "Chromosome size file", true);
		p.addStringArg("-og", "Output file for median enrichment by feature type (requires -ot)", false, null);
		p.addStringArg("-ot", "Output file for table of feature enrichments (requires -og)", false, null);
		p.addStringArg("-om", "Output file for median of max window enrichments per feature type", false, null);
		p.addIntArg("-w", "Window size for max enrichment by feature", false, 50);
		p.addIntArg("-st", "Step size for max window enrichment by feature", false, 10);
		p.addIntArg("-mg", "Max genomic span", false, DEFAULT_MAX_GENOMIC_SPAN);
		p.addBooleanArg("-fc", "Fully contained reads in features", false, DEFAULT_FULLY_CONTAINED);
		p.parse(args);
		String backgroundFile = p.getStringArg("-b");
		String signalFile = p.getStringArg("-s");
		String geneBedFile = p.getStringArg("-g");
		String chrSizeFile = p.getStringArg("-c");
		String outMedian = p.getStringArg("-og");
		String outTable = p.getStringArg("-ot");
		String outMaxWindow = p.getStringArg("-om");
		int windowSize = p.getIntArg("-w");
		int stepSize = p.getIntArg("-st");
		int maxGenomicSpan = p.getIntArg("-mg");
		boolean fullyContained = p.getBooleanArg("-fc");
		
		if(outMedian != null && outTable == null) {
			throw new IllegalArgumentException("If providing -og must also provide -ot.");
		}
		if(outMedian == null && outTable != null) {
			throw new IllegalArgumentException("If providing -ot must also provide -og.");
		}
		
		EnrichmentByFeature ebf = new EnrichmentByFeature(backgroundFile, signalFile, geneBedFile, chrSizeFile, maxGenomicSpan, fullyContained);
		
		if(outMedian != null && outTable != null) ebf.computeAndWriteEnrichments(outMedian, outTable);
		if(outMaxWindow != null) ebf.computeAndWriteMaxWindowEnrichments(outMaxWindow,windowSize,stepSize);
		
		logger.info("All done.");

	}

}
