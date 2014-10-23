package rnaseq;

import guttmanlab.core.util.CommandLineParser;
import guttmanlab.core.util.StringParser;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.TreeSet;

import nextgen.core.annotation.Gene;
import nextgen.core.coordinatesystem.GenomicSpace;
import nextgen.core.coordinatesystem.TranscriptomeSpace;
import nextgen.core.feature.GeneWindow;
import nextgen.core.model.AlignmentModel;
import nextgen.core.model.ScanStatisticDataAlignmentModel;
import nextgen.core.readFilters.FragmentLengthFilter;
import nextgen.core.readFilters.GenomicSpanFilter;
import nextgen.core.readFilters.ProperPairFilter;

import org.apache.log4j.Logger;

import broad.core.math.EmpiricalDistribution;
import broad.pda.annotation.BEDFileParser;

/**
 * @author prussell
 *
 */
public class WindowCountDistribution {

	private ScanStatisticDataAlignmentModel transcriptomeData;
	private ScanStatisticDataAlignmentModel genomeData;
	private boolean hasTranscriptomeData;
	private boolean hasGenomicData;
	private Collection<Gene> genes;
	static Logger logger = Logger.getLogger(WindowCountDistribution.class.getName());
	private TranscriptomeSpace transcriptomeSpace;
	private GenomicSpace genomicSpace;
	private static int DEFAULT_MAX_FRAGMENT_LENGTH = 2000;
	private static int DEFAULT_MAX_GENOMIC_SPAN = 300000;
	private static int DEFAULT_WINDOW_SIZE = 1;
	private static int DEFAULT_STEP_SIZE = 1;
	private static boolean DEFAULT_USE_FRAGMENTS = true;
	private int windowSize;
	private int stepSize;
	
	/**
	 * Instantiate with default parameters
	 * @param bamFile Bam alignments
	 * @param bedFile Bed gene annotation
	 * @throws IOException
	 */
	private WindowCountDistribution(String bamFile, String bedFile, String chrSizeFile) throws IOException {
		this(bamFile, bedFile, chrSizeFile, DEFAULT_MAX_FRAGMENT_LENGTH, DEFAULT_MAX_GENOMIC_SPAN, DEFAULT_WINDOW_SIZE, DEFAULT_STEP_SIZE, DEFAULT_USE_FRAGMENTS);
	}
	
	/**
	 * Instantiate with parameter values
	 * @param bamFile Bam alignments
	 * @param bedFile Bed gene annotation
	 * @param chrSizeFile Chromosome size file
	 * @param maxFragmentLength Max fragment length to count read
	 * @param maxGenomicSpan Max genomic span to count read
	 * @param window Window size for read count distribution
	 * @param step Step size for read count distribution
	 * @throws IOException
	 */
	private WindowCountDistribution(String bamFile, String bedFile, String chrSizeFile, int maxFragmentLength, int maxGenomicSpan, int window, int step, boolean useFragments) throws IOException {

		windowSize = window;
		stepSize = step;
		hasTranscriptomeData = false;
		hasGenomicData = false;
		
		if(bedFile != null) {
			transcriptomeSpace = new TranscriptomeSpace(BEDFileParser.loadDataByChr(new File (bedFile)));
			transcriptomeData = new ScanStatisticDataAlignmentModel(bamFile, transcriptomeSpace, useFragments);
			genes = BEDFileParser.loadData(new File(bedFile));
			// Add read filters
			transcriptomeData.addFilter(new ProperPairFilter());
			transcriptomeData.addFilter(new GenomicSpanFilter(maxGenomicSpan));
			transcriptomeData.addFilter(new FragmentLengthFilter(transcriptomeSpace, maxFragmentLength));
			hasTranscriptomeData = true;
		}
		
		if(chrSizeFile != null) {
			genomicSpace = new GenomicSpace(chrSizeFile);
			genomeData = new ScanStatisticDataAlignmentModel(bamFile, genomicSpace, useFragments);
			// Add read filters
			genomeData.addFilter(new ProperPairFilter());
			genomeData.addFilter(new GenomicSpanFilter(maxGenomicSpan));
			genomeData.addFilter(new FragmentLengthFilter(transcriptomeSpace, maxFragmentLength));
			hasGenomicData = true;
		}
		
		
	}

	/**
	 * Get collection of "RPKM" values for windows within a gene
	 * @param gene The gene
	 * @return Collection of window RPKM values, where the total number of reads for RPKM calculation is reads mapping to the gene
	 */
	private ArrayList<Double> getWindowRpkmValues(Gene gene, AlignmentModel data) {
		double totalCount = data.getCount(gene);
		double rpkmConstant = Math.pow(10, 9) / (windowSize * totalCount);
		logger.info("Gene=" + gene.getName() + "\tWindowSize=" + windowSize + "\tGeneCount=" + totalCount + "\tRpkmConstant=" + rpkmConstant);
		Collection<GeneWindow> windows = gene.getWindows(windowSize, stepSize, 0);
		ArrayList<Double> rtrn = new ArrayList<Double>();
		for(GeneWindow window : windows) {
			rtrn.add(Double.valueOf(rpkmConstant * data.getCount(window)));
		}
		return rtrn;
	}
	
	/**
	 * Get collection of read counts for windows within a gene
	 * @param gene The gene
	 * @return Collection of window read counts
	 */
	private ArrayList<Double> getWindowCounts(Gene gene, AlignmentModel data) {
		Collection<GeneWindow> windows = gene.getWindows(windowSize, stepSize, 0);
		ArrayList<Double> rtrn = new ArrayList<Double>();
		for(GeneWindow window : windows) {
			double count = data.getCount(window);
			rtrn.add(Double.valueOf(count));
			//logger.info(window.getChr() + ":" + window.getStart() + "-" + window.getEnd() + "\t" + count);
			
		}
		return rtrn;
	}
	
	/**
	 * Write raw count for each window of gene to file
	 * @param gene The gene
	 * @param outFile The output file
	 * @throws IOException
	 */
	@SuppressWarnings("unused")
	private void writeWindowCounts(Gene gene, String outFile, AlignmentModel data) throws IOException {
		ArrayList<Double> counts = getWindowCounts(gene, data);
		FileWriter w = new FileWriter(outFile);
		for(Double d : counts) {
			w.write(d.toString() + "\n");
		}
		w.close();
	}

	/**
	 * Write RPKM for each window of gene to file
	 * @param gene The gene
	 * @param outFile The output file
	 * @throws IOException
	 */
	private void writeWindowRpkmValues(Gene gene, String outFile, AlignmentModel data) throws IOException {
		ArrayList<Double> counts = getWindowRpkmValues(gene, data);
		FileWriter w = new FileWriter(outFile);
		for(Double d : counts) {
			w.write(d.toString() + "\n");
		}
		w.close();
	}

	
	/**
	 * Get distribution of "RPKM" values for windows within a gene
	 * @param gene The gene
	 * @param numBins Number of bins for empirical distribution
	 * @return Distribution of window RPKM values, where the total number of reads for RPKM calculation is reads mapping to the gene
	 */
	private EmpiricalDistribution getWindowRpkmDistribution(Gene gene, int numBins, AlignmentModel data) {
		EmpiricalDistribution d = new EmpiricalDistribution(getWindowRpkmValues(gene, data),numBins);
		logger.info("Gene=" + gene.getName() + "\tMedianWindowRpkm=" + d.getMedianOfAllDataValues());
		return d;
	}
	
	/**
	 * Get distribution of read counts for windows within a gene
	 * @param gene The gene
	 * @param numBins Number of bins for empirical distribution
	 * @return Distribution of window read counts
	 */
	private EmpiricalDistribution getWindowCountDistribution(Gene gene, int numBins, AlignmentModel data) {
		EmpiricalDistribution d = new EmpiricalDistribution(getWindowCounts(gene, data),numBins);
		logger.info("Gene=" + gene.getName() + "\tMedianWindowCount=" + d.getMedianOfAllDataValues());
		return d;
	}

	
	/**
	 * Write distribution of read counts for windows within a gene to a file
	 * @param gene The gene
	 * @param numBins Number of bins
	 * @param outFile Output file
	 * @throws IOException
	 */
	@SuppressWarnings("unused")
	private void writeWindowCountDistribution(Gene gene, int numBins, String outFile, AlignmentModel data) throws IOException {
		FileWriter w = new FileWriter(outFile);
		EmpiricalDistribution dist = getWindowCountDistribution(gene, numBins, data);
		w.write("Lower_bound_count\tDensity\n");
		for(int i=0; i < dist.getBinNumber(); i++) {
			w.write(dist.getBinStart(i) + "\t" + dist.getDensity(i) + "\n");
		}
		w.close();
	}
	
	/**
	 * Write distribution of "RPKM" values for windows within a gene to a file
	 * @param gene The gene
	 * @param numBins Number of bins
	 * @param outFile Output file
	 * @throws IOException
	 */
	@SuppressWarnings("unused")
	private void writeWindowRpkmDistribution(Gene gene, int numBins, String outFile, AlignmentModel data) throws IOException {
		FileWriter w = new FileWriter(outFile);
		EmpiricalDistribution dist = getWindowRpkmDistribution(gene, numBins, data);
		w.write("Lower_bound_RPKM\tDensity\n");
		for(int i=0; i < dist.getBinNumber(); i++) {
			w.write(dist.getBinStart(i) + "\t" + dist.getDensity(i) + "\n");
		}
		w.close();
	}

	
	/**
	 * Get the gene set
	 * @return The gene set
	 */
	private Collection<Gene> getGenes() {
		return genes;
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {

		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-b", "Bam file", true);
		p.addStringArg("-g", "Gene bed file for transcriptome space", false, null);
		p.addStringArg("-lg", "List of gene names to write individual files for", false, null);
		p.addStringArg("-cs", "Chromosome size file", false, null);
		p.addStringArg("-lc", "List of entire chromosomes to write individual files for", false, null);
		p.addStringArg("-o", "Output directory", false, ".");
		p.addIntArg("-mf", "Max fragment length to count reads", false, DEFAULT_MAX_FRAGMENT_LENGTH);
		p.addIntArg("-mg", "Max genomic span to count reads", false, DEFAULT_MAX_GENOMIC_SPAN);
		//p.addIntegerArg("-bi", "Number of bins for distribution", false, Integer.valueOf(DEFAULT_NUM_BINS));
		p.addIntArg("-w", "Window size", false, DEFAULT_WINDOW_SIZE);
		p.addIntArg("-s", "Step size", false, DEFAULT_STEP_SIZE);
		p.addBooleanArg("-pe", "Convert paired ends to fragments", false, DEFAULT_USE_FRAGMENTS);
		p.addStringArg("-n", "Sample name", true);
		p.parse(args);
		String geneNameFile = p.getStringArg("-lg");
		String chrSizeFile = p.getStringArg("-cs");
		String chrNameFile = p.getStringArg("-lc");
		String outDir = p.getStringArg("-o");
		int maxFragmentLength = p.getIntArg("-mf");
		int maxGenomicSpan = p.getIntArg("-mg");
		//int numBins = p.getIntegerArg("-bi").intValue();
		String bamFile = p.getStringArg("-b");
		String bedFile = p.getStringArg("-g");
		int windowSize = p.getIntArg("-w");
		boolean fragments = p.getBooleanArg("-pe");
		int stepSize = p.getIntArg("-s");
		String sampleName = p.getStringArg("-n");
		
		// Make the directory
		File dirFile = new File(outDir);
		@SuppressWarnings("unused")
		boolean madeDir = dirFile.mkdir();
		
		WindowCountDistribution pcd = new WindowCountDistribution(bamFile, bedFile, chrSizeFile, maxFragmentLength, maxGenomicSpan, windowSize, stepSize, fragments);
		
		// Find each gene by name and write file
		if(geneNameFile != null) {
			if(!pcd.hasTranscriptomeData) {
				throw new IllegalStateException("Need to provide bed file for transcriptome space.");
			}
			// Read the list of gene names
			FileReader r = new FileReader(geneNameFile);
			BufferedReader b = new BufferedReader(r);
			StringParser s = new StringParser();
			TreeSet<String> genesToCalculate = new TreeSet<String>();
			while(b.ready()) {
				String line = b.readLine();
				s.parse(line);
				if(s.getFieldCount() < 1) continue;
				String geneName = s.asString(0);
				genesToCalculate.add(geneName);
			}
			r.close();
			b.close();
			Collection<Gene> genes = pcd.getGenes();
			for(String geneName : genesToCalculate) {
				boolean found = false;
				logger.info("Calculating distribution for gene " + geneName);
				for(Gene gene : genes) {
					if(!geneName.equals(gene.getName())) continue;
					double geneCount = pcd.transcriptomeData.getCount(gene);
					logger.info("Total gene count = " + geneCount);
					String outRpkm = outDir + "/window_rpkm_" + sampleName + "_" + geneName;
					pcd.writeWindowRpkmValues(gene, outRpkm, pcd.transcriptomeData);
					found = true;
					continue;
				}
				if(!found) {
					logger.warn("Gene " + geneName + " does not exist in transcriptome space.");
				}
			}
		}
		
		// Find each chromosome by name and write file
		if(chrNameFile != null) {
			if(!pcd.hasGenomicData) {
				throw new IllegalStateException("Need to provide chromosome size file for genome space.");
			}
			// Read the list of chromosome names
			FileReader r = new FileReader(chrNameFile);
			BufferedReader b = new BufferedReader(r);
			StringParser s = new StringParser();
			TreeSet<String> chrsToCalculate = new TreeSet<String>();
			while(b.ready()) {
				String line = b.readLine();
				s.parse(line);
				if(s.getFieldCount() < 1) continue;
				String chrName = s.asString(0);
				chrsToCalculate.add(chrName);
			}
			r.close();
			b.close();
			for(String chrName : chrsToCalculate) {
				logger.info("Calculating distribution for chromosome " + chrName);
				Gene chrAsGene = new Gene(pcd.genomicSpace.getReferenceAnnotation(chrName));
				double chrCount = pcd.genomeData.getCount(chrAsGene);
				logger.info("Total chr count = " + chrCount);
				String outRpkm = outDir + "/window_rpkm_" + sampleName + "_" + chrName;
				pcd.writeWindowRpkmValues(chrAsGene, outRpkm, pcd.genomeData);
				continue;
			}
		}

		
		
	}

}
