/**
 * 
 */
package annotation;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.log4j.Logger;

import broad.core.math.EmpiricalDistribution;
import broad.core.math.Statistics;
import broad.core.parser.CommandLineParser;
import broad.core.parser.StringParser;
import broad.pda.annotation.BEDFileParser;
import broad.pda.seq.protection.SampleData;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.BasicAnnotation;
import nextgen.core.annotation.Gene;
import nextgen.core.utils.AnnotationUtils;

/**
 * @author prussell
 *
 */
public class FeatureSetStatistics {

	private Map<String, Collection<Gene>> genesByChr;
	private Map<String, TreeSet<Gene>> featuresByChr;
	protected Map<String, Gene> genesByName;
	private Map<String, Gene> featuresByName;
	private Map<Gene, Collection<Gene>> featuresByGene;
	private String outputDir;
	private boolean useExistingFiles;
	protected String bedAnnotation;

	protected Map<ExpressionQuantileBin, EmpiricalDistribution> distNumFeaturesPerGene;
	protected Map<ExpressionQuantileBin, EmpiricalDistribution> distPctFeaturesPerGene;
	protected EmpiricalDistribution distFeatureSize;
	protected Map<ExpressionQuantileBin, EmpiricalDistribution> distAvgDistBetweenFeatures;
	protected Map<ExpressionQuantileBin, EmpiricalDistribution> distMedianFeatureSizeByGene;
	protected Map<ExpressionQuantileBin, EmpiricalDistribution> distMaxDistanceBtwConsecutiveFeatures;
	protected Map<ExpressionQuantileBin, EmpiricalDistribution> distMedianDistanceBtwConsecutiveFeatures;
	
	private SampleData expressionData;
	protected double[] sortedExpressionLevels;
	protected Map<Gene, Double> expressionLevels;
	private Map<Gene, Double> expressionQuantiles;
	private boolean checkExpression;
	private ArrayList<ExpressionQuantileBin> expressionQuantileBins;
	private static int DEFAULT_NUM_QUANTILE_BINS = 4;
	
	private StatisticAverageDistanceBetweenFeatures statAvgDistBtwFeatures;
	private StatisticFeatureSize statFeatureSize;
	private StatisticMaxDistanceBetweenConsecutiveFeatures statMaxDistBtwConsecutiveFeatures;
	private StatisticMedianDistanceBetweenConsecutiveFeatures statMedianDistBtwConsecutiveFeatures;
	private StatisticNumFeaturesPerGene statNumFeaturesPerGene;
	private StatisticPercentFeatureCoveragePerGene statPctFeaturesPerGene;
	private StatisticMedianFeatureSizeByGene statMedianFeatureSizeByGene;
	
	private static boolean DEFAULT_FIRST_READ_TRANSCRIPTION_STRAND = false;
	
	protected static Logger logger = Logger.getLogger(FeatureSetStatistics.class.getName());
	
	private FeatureSetStatistics(String bedGenes, String bedFeatures, String outDir, boolean readFromExistingFiles) throws IOException {
		this(bedGenes, bedFeatures, outDir, readFromExistingFiles, null, 1, DEFAULT_NUM_QUANTILE_BINS, DEFAULT_FIRST_READ_TRANSCRIPTION_STRAND);
	}
	
	private FeatureSetStatistics(String bedGenes, String bedFeatures, String outDir, boolean readFromExistingFiles, String expressionBamFile, double expressionPvalCutoff, int numQuantileBins, boolean firstReadTranscriptionStrand) throws IOException {
		
		bedAnnotation = bedGenes;
		
		outputDir = outDir;
		File outDirFile = new File(outputDir);
		@SuppressWarnings("unused")
		boolean madeDir = outDirFile.mkdir();
		if(!outDirFile.exists()) {
			throw new IOException("Could not create directory " + outputDir);
		}
		
		useExistingFiles = readFromExistingFiles;
		
		genesByChr = new TreeMap<String, Collection<Gene>>();
		genesByName = new TreeMap<String, Gene>();
		Map<String, Collection<Gene>> genesFromBed = BEDFileParser.loadDataByChr(new File(bedAnnotation));
		for(String chr : genesFromBed.keySet()) {
			TreeSet<Gene> genesThisChr = new TreeSet<Gene>();
			genesThisChr.addAll(genesFromBed.get(chr));
			genesByChr.put(chr, genesThisChr);
			for(Gene gene : genesThisChr) {
				genesByName.put(gene.getName(), gene);
			}
		}
		
		featuresByChr = new TreeMap<String, TreeSet<Gene>>();
		featuresByName = new TreeMap<String, Gene>();
		Map<String, Collection<Gene>> featuresFromBed = BEDFileParser.loadDataByChr(new File(bedFeatures));
		for(String chr : featuresFromBed.keySet()) {
			TreeSet<Gene> featuresThisChr = new TreeSet<Gene>();
			featuresThisChr.addAll(featuresFromBed.get(chr));
			featuresByChr.put(chr, featuresThisChr);
			for(Gene feature : featuresThisChr) {
				featuresByName.put(feature.getName(), feature);
			}
		}
		
		logger.info("Done loading genes and features.");
		
		featuresByGene = new TreeMap<Gene, Collection<Gene>>();
		String cachedFeatureGeneMappingFile = getCachedFeatureGeneMappingFileName(bedAnnotation, bedFeatures);
		mapFeaturesToGenes(cachedFeatureGeneMappingFile);
		
		expressionQuantileBins = new ArrayList<ExpressionQuantileBin>();
		expressionQuantileBins.add(allQuantiles());
		expressionQuantileBins.addAll(getQuantileBinsByNumBins(numQuantileBins));
		
		distNumFeaturesPerGene = new HashMap<ExpressionQuantileBin, EmpiricalDistribution>();
		distPctFeaturesPerGene = new HashMap<ExpressionQuantileBin, EmpiricalDistribution>();
		distFeatureSize = null;
		distAvgDistBetweenFeatures = new HashMap<ExpressionQuantileBin, EmpiricalDistribution>();
		distMedianFeatureSizeByGene = new HashMap<ExpressionQuantileBin, EmpiricalDistribution>();
		distMaxDistanceBtwConsecutiveFeatures = new HashMap<ExpressionQuantileBin, EmpiricalDistribution>();
		distMedianDistanceBtwConsecutiveFeatures = new HashMap<ExpressionQuantileBin, EmpiricalDistribution>();
		
		statAvgDistBtwFeatures = new StatisticAverageDistanceBetweenFeatures();
		statFeatureSize = new StatisticFeatureSize();
		statMaxDistBtwConsecutiveFeatures = new StatisticMaxDistanceBetweenConsecutiveFeatures();
		statMedianDistBtwConsecutiveFeatures = new StatisticMedianDistanceBetweenConsecutiveFeatures();
		statNumFeaturesPerGene = new StatisticNumFeaturesPerGene();
		statPctFeaturesPerGene = new StatisticPercentFeatureCoveragePerGene();
		statMedianFeatureSizeByGene = new StatisticMedianFeatureSizeByGene();
		
		if(expressionBamFile == null) {
			checkExpression = false;
			logger.warn("Not checking gene expression; counting all genes.");
		} else {
			checkExpression = true;
			logger.info("Loading expression data...");
			expressionData = new SampleData(expressionBamFile, firstReadTranscriptionStrand, genesByChr, 0, 0, expressionPvalCutoff, true, true);
			computeExpressionLevels();
			computeExpressionQuantiles();
			logger.info("Done loading expression data.");
		}
		
	}
	
	/**
	 * Compute and store expression levels of genes
	 * @throws IOException 
	 */
	private void computeExpressionLevels() throws IOException {
		
		logger.info("");
		logger.info("Computing expression levels...");
		
		ExpressionScoreFile e = new ExpressionScoreFile();
		if(e.readFromFile()) {
			return;
		}
		
		ArrayList<Double> levelsList = new ArrayList<Double>();
		expressionLevels = new TreeMap<Gene, Double>();
		int numDone = 0;
		for(String name : genesByName.keySet()) {
			Gene gene = genesByName.get(name);
			double p = expressionData.getGeneScanPval(gene);
			double level = -1 * Math.log(p);
			expressionLevels.put(gene, Double.valueOf(level));
			levelsList.add(Double.valueOf(level));
			numDone++;
			if(numDone % 1000 == 0) {
				logger.info("Finished " + numDone + " genes.");
			}
		}
		sortedExpressionLevels = new double[levelsList.size()];
		for(int i=0; i < sortedExpressionLevels.length; i++) {
			sortedExpressionLevels[i] = levelsList.get(i).doubleValue();
		}
		Arrays.sort(sortedExpressionLevels);
		
		e.writeToFile();
		
	}
	
	/**
	 * Compute and store expression quantiles of genes
	 */
	private void computeExpressionQuantiles() {
		logger.info("Computing expression quantiles...");
		expressionQuantiles = new TreeMap<Gene, Double>();
		for(String name : genesByName.keySet()) {
			Gene gene = genesByName.get(name);
			//TODO bug might be here
			double quantile = Math.abs((Arrays.binarySearch(sortedExpressionLevels, expressionLevels.get(gene).doubleValue()))) / (double) sortedExpressionLevels.length;
			expressionQuantiles.put(gene, Double.valueOf(quantile));
		}
		logger.info("Done computing expression quantiles.");
	}
	
	/**
	 * Whether the gene is expressed
	 * @param gene The gene
	 * @return True if gene is expressed or not checking expression
	 */
	private boolean isExpressed(Gene gene) {
		if(!checkExpression) return true;
		return expressionData.isExpressed(gene);
	}
	
	
	
	/**
	 * Identify each gene with set of overlapping features
	 * @throws IOException 
	 */
	private void mapFeaturesToGenes(String cachedFeatureGeneMappingFile) throws IOException {
		
		logger.info("");
		logger.info("Identifying features with genes...");
		
		if(useExistingFiles) {
			if(readFeatureGeneMappingFromCachedFile(cachedFeatureGeneMappingFile)) {
				return;
			}
		}
		
		featuresByGene.clear();
		for(String chr : genesByChr.keySet()) {
			int numGenes = 0;
			int numFeatures = 0;
			for(Gene gene : genesByChr.get(chr)) {
				Collection<Gene> geneFeatures = new TreeSet<Gene>();
				if (featuresByChr.containsKey(chr)) {
					for(Gene feature : featuresByChr.get(chr)) {
						if(feature.overlaps(gene)) {
							geneFeatures.add(feature);
						}
					}
				}				
				featuresByGene.put(gene, geneFeatures);
				int f = featuresByGene.get(gene).size();
				if(f > 0) {
					numGenes++;
					numFeatures += f;
				}
			}
			logger.info(chr + " " + numGenes + "/" + genesByChr.get(chr).size() + " genes have " + numFeatures + " features.");
		}
		
		logger.info("Done identifying features with genes.");
		
		if(useExistingFiles) {
			writeFeatureGeneMapping(cachedFeatureGeneMappingFile);
		}
		
	}
	
	/**
	 * Construct name of cached file for mapping of features to genes
	 * @param bedAnnotation Annotation bed file
	 * @param bedFeatures Feature bed file
	 * @return Name of cached file
	 */
	private static String getCachedFeatureGeneMappingFileName(String bedAnnotation, String bedFeatures) {
		
		String rtrn = "gene_feature_mapping_";
		
		if(bedAnnotation.endsWith(".bed")) {
			rtrn += bedAnnotation.substring(0, bedAnnotation.length() - 4) + "_";
		} else {
			rtrn += bedAnnotation + "_";
		}
		
		if(bedFeatures.endsWith(".bed")) {
			rtrn += bedFeatures.substring(0, bedFeatures.length() - 4);
		} else {
			rtrn += bedFeatures;
		}
		
		return rtrn;
		
	}
	
	/**
	 * Write feature to gene mapping to file
	 * @param cachedFeatureGeneMappingFile
	 * @throws IOException
	 */
	private void writeFeatureGeneMapping(String cachedFeatureGeneMappingFile) throws IOException {
		
		logger.info("");
		logger.info("Writing feature to gene mapping to file " + cachedFeatureGeneMappingFile + "...");
		
		FileWriter w = new FileWriter(cachedFeatureGeneMappingFile);
		for(Gene gene : featuresByGene.keySet()) {
			String geneName = gene.getName();
			for(Gene feature : featuresByGene.get(gene)) {
				String featureName = feature.getName();
				w.write(geneName + "\t" + featureName + "\n");
			}
		}
		w.close();
		
		logger.info("Done writing to file.");
		
	}
	
	/**
	 * Read feature to gene mapping from file and store
	 * @param cachedFeatureGeneMappingFile Cached file
	 * @return Whether the operation was successful
	 */
	private boolean readFeatureGeneMappingFromCachedFile(String cachedFeatureGeneMappingFile) {
		
		File cachedFile = new File(cachedFeatureGeneMappingFile);
		if(!cachedFile.exists()) return false;
		
		logger.info("Trying to read feature to gene mapping from file " + cachedFeatureGeneMappingFile + "...");
		
		try {
			int numRead = 0;
			featuresByGene.clear();
			for(String chr : genesByChr.keySet()) {
				for(Gene gene : genesByChr.get(chr)) {
					Collection<Gene> geneFeatures = new TreeSet<Gene>();
					featuresByGene.put(gene, geneFeatures);
				}
			}
			FileReader r = new FileReader(cachedFeatureGeneMappingFile);
			BufferedReader b = new BufferedReader(r);
			StringParser s = new StringParser();
			while(b.ready()) {
				String line = b.readLine();
				numRead++;
				if(numRead % 1000 == 0) {
					logger.info("Read " + numRead + " mappings.");
				}
				s.parse(line);
				String geneName = s.asString(0);
				String featureName = s.asString(1);
				if(!genesByName.containsKey(geneName)) {
					logger.warn("Gene name " + geneName + " not recognized. Skipping mapping: " + line);
					continue;
				}
				if(!featuresByName.containsKey(featureName)) {
					logger.warn("Feature name " + featureName + " not recognized. Skipping mapping: " + line);
					continue;
				}
				Gene gene = genesByName.get(geneName);
				Gene feature = featuresByName.get(featureName);
				if(!featuresByGene.containsKey(gene)) {
					Collection<Gene> geneFeatures = new TreeSet<Gene>();
					featuresByGene.put(gene, geneFeatures);
				}
				featuresByGene.get(gene).add(feature);
			}
			r.close();
			b.close();
			
		} catch(Exception e) {
			logger.info("Could not read gene to feature mapping from file: " + cachedFeatureGeneMappingFile);
			logger.info("Due to exception:");
			e.printStackTrace();
			return false;
		}
		
		logger.info("Done reading from file.");
		
		return true;
		
	}
	
	/**
	 * Write distributions of number of features per gene to files in output directory
	 * @throws IOException
	 */
	private void writeDistributions(Statistic stat) throws IOException {
		String statName = stat.toString();
		logger.info("");
		logger.info("Writing distributions of " + statName + "...");
		Map<ExpressionQuantileBin, String> outFiles = new HashMap<ExpressionQuantileBin, String>();
		for(ExpressionQuantileBin bin : expressionQuantileBins) {
			String outFile = outputDir + "/dist_" + statName + "_" + bin.toString();
			outFiles.put(bin, outFile);
		}
		boolean allOutFilesExist = true;
		for(ExpressionQuantileBin bin : outFiles.keySet()) {
			File of = new File(outFiles.get(bin));
			if(!of.exists()) {
				allOutFilesExist = false;
				break;
			}
		}
		if (allOutFilesExist && useExistingFiles) {
			logger.info("All output files already exist. Not recomputing distributions.");
			return;
		}		
		if (stat.getDistributions() == null) {
			stat.calculate();
		}
		if (stat.getDistributions().isEmpty()) {
			stat.calculate();
		}
		for (ExpressionQuantileBin bin : expressionQuantileBins) {
			EmpiricalDistribution dist = stat.getDistributions().get(bin);
			if(dist == null) {
				logger.warn("Statistic " + stat.toString() + " does not have distribution for expression in bin " + bin.toString());
				continue;
			}
			String outFile = outFiles.get(bin);
			logger.info("Writing distribution for genes with expression in bin " + bin.toString() + " to file " + outFile);
			dist.write(outFile);
		}
		logger.info("Done writing files.");
	}
	
	/**
	 * Make distribution of number of features per gene
	 * @throws IOException 
	 */
	protected void makeFeatureCoverageDistributions() throws IOException {
		
		logger.info("");
		logger.info("Making distributions of feature coverage statistics...");
		
		// Number of features per gene
		String outFileNumFeaturesPerGene = outputDir + "/table_num_features_per_gene";
		FileWriter wn = new FileWriter(outFileNumFeaturesPerGene);
		String headerN = "gene_name\tgene_size\texpression_quantile\tnum_features";
		wn.write(headerN + "\n");

		Map<ExpressionQuantileBin, Collection<Double>> countPerGene = new HashMap<ExpressionQuantileBin, Collection<Double>>();
		for(ExpressionQuantileBin bin : expressionQuantileBins) {
			countPerGene.put(bin, new ArrayList<Double>());
		}

		// Percentage of each gene covered by features
		String outFilePctFeaturesPerGene = outputDir + "/table_pct_features_per_gene";
		FileWriter wp = new FileWriter(outFilePctFeaturesPerGene);
		String headerP = "gene_name\tgene_size\texpression_quantile\tpct_covered";
		wp.write(headerP + "\n");

		Map<ExpressionQuantileBin, Collection<Double>> pctPerGene = new HashMap<ExpressionQuantileBin, Collection<Double>>();
		for(ExpressionQuantileBin bin : expressionQuantileBins) {
			pctPerGene.put(bin, new ArrayList<Double>());
		}
		
		int numNotExpressed = 0;
		for(String chr : genesByChr.keySet()) {
			logger.info(chr);
			for(Gene gene : genesByChr.get(chr)) {
				
				try {
					// Check gene expression
					if(!isExpressed(gene)) {
						numNotExpressed++;
						continue;
					}
					
					double quantile = expressionQuantiles.get(gene).doubleValue();
					
					// Make stats
					int count = featuresByGene.get(gene).size();
					int geneSize = gene.getSize();
					double p = (double)getNumBasesCoveredByFeatures(gene) / (double)geneSize;
					
					for (ExpressionQuantileBin bin : expressionQuantileBins) {
						if(bin.inBin(quantile)) {
							// Add stats
							pctPerGene.get(bin).add(Double.valueOf(p));
							countPerGene.get(bin).add(Double.valueOf(count));
						}
					}
					// Write to tables
					wn.write(gene.getName() + "\t" + geneSize + "\t" + quantile + "\t" + count + "\n");
					wp.write(gene.getName() + "\t" + geneSize + "\t" + quantile + "\t" + p + "\n");
				} catch (NullPointerException e) {
					logger.warn("Caught exception when getting stats for gene " + gene.getName());
					e.printStackTrace();
				}
				
			}
		}
		
		logger.warn("Skipped " + numNotExpressed + " genes that are not expressed.");
		
		for(ExpressionQuantileBin bin : expressionQuantileBins) {
			int numBins = Math.max((int) Statistics.max(countPerGene.get(bin)),1);
			distNumFeaturesPerGene.put(bin, new EmpiricalDistribution(countPerGene.get(bin), numBins));
			distPctFeaturesPerGene.put(bin, new EmpiricalDistribution(pctPerGene.get(bin), 50));
		}
		wn.close();
		wp.close();
		logger.info("Done making distributions of feature coverage statistics.");
		
	}
	
	/**
	 * Make distributions of feature distribution within gene
	 * @throws IOException 
	 */
	protected void makeFeatureStatsDistributions() throws IOException {
		
		logger.info("");
		logger.info("Making distributions of feature stats by gene...");
		
		// Median feature size by gene
		String outFileMedianSize = outputDir + "/table_median_feature_size_by_gene";
		FileWriter wm = new FileWriter(outFileMedianSize);
		String headerM = "gene_name\tnum_features\tmedian_feature_size";
		wm.write(headerM + "\n");
		
		// Average distance between features
		String outFileAvgDist = outputDir + "/table_avg_dist_between_features";
		FileWriter wa = new FileWriter(outFileAvgDist);
		String headerA = "gene_name\tgene_size\tnum_features\tavg_dist_between_features";
		wa.write(headerA + "\n");

		Map<ExpressionQuantileBin, Collection<Double>> avgDistances = new HashMap<ExpressionQuantileBin, Collection<Double>>();
		for(ExpressionQuantileBin bin : expressionQuantileBins) {
			avgDistances.put(bin, new ArrayList<Double>());
		}
		
		
		Map<ExpressionQuantileBin, Collection<Double>> medianSizePerGene = new HashMap<ExpressionQuantileBin, Collection<Double>>();
		for(ExpressionQuantileBin bin : expressionQuantileBins) {
			medianSizePerGene.put(bin, new ArrayList<Double>());
		}
		
		int numNotExpressed = 0;
		int numNoFeatures = 0;
		for(String chr : genesByChr.keySet()) {
			logger.info(chr);
			for(Gene gene : genesByChr.get(chr)) {
				
				try {
					// Check gene expression
					if(!isExpressed(gene)) {
						numNotExpressed++;
						continue;
					}
					
					double quantile = expressionQuantiles.get(gene).doubleValue();
					
					// Check if gene has features
					Collection<Gene> features = featuresByGene.get(gene);
					int numFeatures = features.size();
					if(numFeatures == 0) {
						numNoFeatures++;
						continue;
					}
					
					// Make stats
					int len = gene.getSize();
					int numBasesInFeatures = getNumBasesCoveredByFeatures(gene);
					int basesNotCovered = len - numBasesInFeatures;
					double avgDistBtwFeatures = (double)basesNotCovered / (double)numFeatures;
					Collection<Double> sizes = new ArrayList<Double>();
					for(Gene feature : features) {
						sizes.add(Double.valueOf(feature.getSize()));
					}
					double medianSize = Statistics.median(sizes);
					
					// Add stats
					for (ExpressionQuantileBin bin : expressionQuantileBins) {
						if(bin.inBin(quantile)) {
							// Add stats
							avgDistances.get(bin).add(Double.valueOf(avgDistBtwFeatures));
							medianSizePerGene.get(bin).add(Double.valueOf(medianSize));
						}
					}
					
					// Write to tables
					wa.write(gene.getName() + "\t" + gene.getSize() + "\t" + numFeatures + "\t" + avgDistBtwFeatures + "\n");
					wm.write(gene.getName() + "\t" + numFeatures + "\t" + medianSize + "\n");
				} catch (NullPointerException e) {
					logger.warn("Caught exception when getting stats for gene " + gene.getName());
					e.printStackTrace();
				}
			}
		}
		logger.warn("Skipped " + numNotExpressed + " genes that are not expressed.");
		logger.warn("Skipped " + numNoFeatures + " genes that have no features.");
		for(ExpressionQuantileBin bin : expressionQuantileBins) {
			distMedianFeatureSizeByGene.put(bin, new EmpiricalDistribution(medianSizePerGene.get(bin), 100));
			distAvgDistBetweenFeatures.put(bin, new EmpiricalDistribution(avgDistances.get(bin), 50));
		}
		
		wa.close();
		wm.close();
		
		logger.info("Done making distributions of feature stats per gene.");		
		
	}

	/**
	 * Make distributions of consecutive features within gene
	 * @throws IOException 
	 */
	protected void makeConsecutiveFeatureDistributions() throws IOException {
		
		logger.info("");
		logger.info("Making distributions of consecutive feature stats by gene...");
		
		// Median distance between consecutive features
		String outFileMedianDist = outputDir + "/table_median_dist_between_consecutive_features_by_gene";
		FileWriter wm = new FileWriter(outFileMedianDist);
		String headerM = "gene_name\tnum_features\tmedian_dist_between_consecutive_features";
		wm.write(headerM + "\n");
		
		// Longest distance between consecutive featuers
		String outFileLongestDist = outputDir + "/table_max_dist_between_consecutive_features_by_gene";
		FileWriter wl = new FileWriter(outFileLongestDist);
		String headerL = "gene_name\tnum_features\tlongest_dist_between_consecutive_features";
		wl.write(headerL + "\n");

		
		Map<ExpressionQuantileBin, Collection<Double>> medianDistances = new HashMap<ExpressionQuantileBin, Collection<Double>>();
		for(ExpressionQuantileBin bin : expressionQuantileBins) {
			medianDistances.put(bin, new ArrayList<Double>());
		}
		
		
		Map<ExpressionQuantileBin, Collection<Double>> longestDistances = new HashMap<ExpressionQuantileBin, Collection<Double>>();
		for(ExpressionQuantileBin bin : expressionQuantileBins) {
			longestDistances.put(bin, new ArrayList<Double>());
		}

		
		int numNotExpressed = 0;
		int numLessThanTwoFeatures = 0;
		for(String chr : genesByChr.keySet()) {
			logger.info(chr);
			for(Gene gene : genesByChr.get(chr)) {
				
				try {
					// Check gene expression
					if(!isExpressed(gene)) {
						numNotExpressed++;
						continue;
					}
					
					double quantile = expressionQuantiles.get(gene).doubleValue();
					
					// Check if gene has features
					Collection<Gene> features = featuresByGene.get(gene);
					int numFeatures = features.size();
					if(numFeatures < 2) {
						numLessThanTwoFeatures++;
						continue;
					}
					
					// Make stats
					Collection<Double> distancesBtwConsecutiveFeatures = new ArrayList<Double>();
					ArrayList<Gene> featuresInOrder = new ArrayList<Gene>();
					for(Gene feature : features) {
						featuresInOrder.add(feature);
					}
					for(int i=0; i<featuresInOrder.size() - 1; i++) {
						Gene feature1 = featuresInOrder.get(i);
						Gene feature2 = featuresInOrder.get(i+1);
						double dist = transcriptDistanceBetweenFeatures(gene, feature1, feature2);
						distancesBtwConsecutiveFeatures.add(Double.valueOf(dist));
					}
					double median = Statistics.median(distancesBtwConsecutiveFeatures);
					double max = Statistics.max(distancesBtwConsecutiveFeatures);
					
					// Add stats
					for (ExpressionQuantileBin bin : expressionQuantileBins) {
						if(bin.inBin(quantile)) {
							// Add stats
							medianDistances.get(bin).add(Double.valueOf(median));
							longestDistances.get(bin).add(Double.valueOf(max));
						}
					}

					
					// Write to tables
					wl.write(gene.getName() + "\t" + numFeatures + "\t" + max + "\n");
					wm.write(gene.getName() + "\t" + numFeatures + "\t" + median + "\n");
				} catch (NullPointerException e) {
					logger.warn("Caught exception when getting stats for gene " + gene.getName());
					e.printStackTrace();
				}
			}
		}
		logger.warn("Skipped " + numNotExpressed + " genes that are not expressed.");
		logger.warn("Skipped " + numLessThanTwoFeatures + " genes that have fewer than 2 features.");
		for(ExpressionQuantileBin bin : expressionQuantileBins) {
			distMedianDistanceBtwConsecutiveFeatures.put(bin, new EmpiricalDistribution(medianDistances.get(bin), 50));
			distMaxDistanceBtwConsecutiveFeatures.put(bin, new EmpiricalDistribution(longestDistances.get(bin), 50));
		}
		
		wl.close();
		wm.close();
		
		logger.info("Done making distributions of consecutive feature stats per gene.");
		
	}

	
	/**
	 * Get distance along transcript between two features
	 * @param parentTranscript Parent transcript
	 * @param feature1 One feature in no particular order
	 * @param feature2 Another feature in no particular order
	 * @return The distance between the end of the left feature and the start of the right feature, or zero if feature spans overlap
	 */
	private static int transcriptDistanceBetweenFeatures(Gene parentTranscript, Gene feature1, Gene feature2) {
		
		if(!parentTranscript.overlaps(feature1)) {
			throw new IllegalArgumentException("Gene " + parentTranscript.getName() + " does not overlap feature " + feature1.getName());
		}
		if(!parentTranscript.overlaps(feature2)) {
			throw new IllegalArgumentException("Gene " + parentTranscript.getName() + " does not overlap feature " + feature2.getName());
		}
		
		Annotation span1 = new BasicAnnotation(feature1.getChr(), feature1.getStart(), feature1.getEnd());
		Annotation span2 = new BasicAnnotation(feature2.getChr(), feature2.getStart(), feature2.getEnd());
		Annotation intersect1 = parentTranscript.intersect(span1);
		Annotation intersect2 = parentTranscript.intersect(span2);
		if(intersect1.overlaps(intersect2)) {
			return 0;
		}
		
		Gene leftGene = feature1.compareTo(feature2) < 0 ? feature1 : feature2;
		Gene rightGene = feature1.compareTo(feature2) < 0 ? feature2 : feature1;
		int leftGeneEnd = leftGene.getEnd();
		int rightGeneStart = rightGene.getStart();
		int leftGeneEndTranscript = parentTranscript.genomicToTranscriptPosition(leftGeneEnd);
		int rightGeneStartTranscript = parentTranscript.genomicToTranscriptPosition(rightGeneStart);
		
		int rtrn = Math.abs(rightGeneStartTranscript - leftGeneEndTranscript);
		return rtrn;
	}
	
	/**
	 * Get number of bases in gene covered by features
	 * Bases covered by more than one feature are counted once
	 * @param gene The gene
	 * @return Number of bases covered by features
	 */
	private int getNumBasesCoveredByFeatures(Gene gene) {
		Collection<Gene> features = featuresByGene.get(gene);
		TreeSet<Annotation> featureTree = new TreeSet<Annotation>();
		featureTree.addAll(features);
		Collection<Annotation> mergedFeatures = AnnotationUtils.mergeOverlappingBlocks(featureTree);
		int numOverlappingBases = 0;
		for(Annotation feature : mergedFeatures) {
			numOverlappingBases += gene.getOverlap(feature);
		}
		return numOverlappingBases;
	}
	
	/**
	 * Make distribution of feature sizes
	 * @throws IOException 
	 */
	protected void makeDistFeatureSize() throws IOException {
		
		logger.info("");
		logger.info("Making distribution of feature sizes...");
		
		String outFile = outputDir + "/table_feature_length";
		FileWriter w = new FileWriter(outFile);
		String header = "feature_name\tfeature_length";
		w.write(header + "\n");
		
		Collection<Double> sizes = new ArrayList<Double>();
		for(String chr : featuresByChr.keySet()) {
			logger.info(chr);
			for(Gene feature : featuresByChr.get(chr)) {
				int len = feature.size();
				w.write(feature.getName() + "\t" + len + "\n");
				sizes.add(Double.valueOf(len));
			}
		}
		distFeatureSize = new EmpiricalDistribution(sizes, 5000);
		
		w.close();
		
		logger.info("Done making distribution of feature lengths.");		
		
	}
	
	private void writeSummary() throws IOException {
		
		String outFile = outputDir + "/summary";
		
		logger.info("");
		logger.info("Writing summary to file " + outFile + "...");
		FileWriter w = new FileWriter(outFile);
		int numExpressed = 0;
		int numWithFeatures = 0;
		for(String name : genesByName.keySet()) {
			Gene gene = genesByName.get(name);
			if(isExpressed(gene)) {
				numExpressed++;
			}
			if(featuresByGene.get(gene).size() > 0) {
				numWithFeatures++;
			}
		}
		w.write("Number of genes: " + genesByName.size() + "\n");
		w.write("Number of expressed genes: " + numExpressed + "\n");
		w.write("Number of genes with features: " + numWithFeatures + "\n");
		w.write("Number of features: " + featuresByName.size() + "\n");
		w.write("\n");
		String header = "distribution\tunit\texpression_quantile_bin\tmedian";
		w.write(header + "\n");
		w.write("feature_size\tall_features\t" + allQuantiles().toString() + "\t" + distFeatureSize.getMedianOfAllDataValues() + "\n");
		for(ExpressionQuantileBin bin : expressionQuantileBins) {
			w.write("num_features_per_gene\tall_genes\t" + bin.toString() + "\t" + distNumFeaturesPerGene.get(bin).getMedianOfAllDataValues() + "\n");
			w.write("percent_of_gene_covered_by_features\tall_genes\t" + bin.toString() + "\t" + distPctFeaturesPerGene.get(bin).getMedianOfAllDataValues() + "\n");
			w.write("median_feature_size_by_gene\tgenes_with_at_least_one_feature\t" + bin.toString() + "\t" + distMedianFeatureSizeByGene.get(bin).getMedianOfAllDataValues() + "\n");
			w.write("average_distance_between_features\tgenes_with_at_least_one_feature\t" + bin.toString() + "\t" + distAvgDistBetweenFeatures.get(bin).getMedianOfAllDataValues() + "\n");
			w.write("max_distance_between_consecutive_features\tgenes_with_at_least_two_features\t" + bin.toString() + "\t" + distMaxDistanceBtwConsecutiveFeatures.get(bin).getMedianOfAllDataValues() + "\n");
			w.write("median_distance_between_consecutive_features\tgenes_with_at_least_two_features\t" + bin.toString() + "\t" + distMedianDistanceBtwConsecutiveFeatures.get(bin).getMedianOfAllDataValues() + "\n");
		}
		logger.info("Done writing summary.");
		
		w.close();
	}
	
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-g", "Gene bed file", true);
		p.addStringArg("-f", "Feature bed file", true);
		p.addStringArg("-o", "Output directory", true);
		p.addBooleanArg("-e", "Try to use existing files from previous run", false, false);
		p.addStringArg("-eb", "Bam file for gene expression", false, null);
		p.addDoubleArg("-ep", "P-value cutoff for gene expression (required to check gene expression)", false, 0.01);
		p.addIntArg("-nb", "Number of expression quantile bins (>=1)", false, DEFAULT_NUM_QUANTILE_BINS);
		p.addBooleanArg("-ft", "First read is transcription strand", false, DEFAULT_FIRST_READ_TRANSCRIPTION_STRAND);
		p.parse(args);
		String geneBed = p.getStringArg("-g");
		String featureBed = p.getStringArg("-f");
		String outDir = p.getStringArg("-o");
		boolean useExisting = p.getBooleanArg("-e");
		String expressionBam = p.getStringArg("-eb");
		double expressionPvalCutoff = p.getDoubleArg("-ep");
		int numQuantileBins = p.getIntArg("-nb");
		boolean firstReadTranscriptionStrand = p.getBooleanArg("-ft");
		
		FeatureSetStatistics f = new FeatureSetStatistics(geneBed, featureBed, outDir, useExisting, expressionBam, expressionPvalCutoff, numQuantileBins, firstReadTranscriptionStrand);
		
		f.writeDistributions(f.statAvgDistBtwFeatures);
		f.writeDistributions(f.statFeatureSize);
		f.writeDistributions(f.statMaxDistBtwConsecutiveFeatures);
		f.writeDistributions(f.statMedianDistBtwConsecutiveFeatures);
		f.writeDistributions(f.statNumFeaturesPerGene);
		f.writeDistributions(f.statPctFeaturesPerGene);
		f.writeDistributions(f.statMedianFeatureSizeByGene);
		f.writeSummary();
		
		logger.info("");
		logger.info("All done.");
		
	}
	
	/**
	 * Get list of all bins of specified size
	 * @param binSize Bin size
	 * @return List of bins of this size
	 */
	public ArrayList<ExpressionQuantileBin> getQuantileBinsBySize(double binSize) {
		if(binSize <= 0 || binSize > 1) {
			throw new IllegalArgumentException("Bin size must be between 0 and 1");
		}
		ArrayList<ExpressionQuantileBin> rtrn = new ArrayList<ExpressionQuantileBin>();
		double min = 0;
		double max = binSize;
		while(max <= 1) {
			rtrn.add(new ExpressionQuantileBin(min,max));
			min = max;
			max = min + binSize;
		}
		if(min < 1) {
			rtrn.add(new ExpressionQuantileBin(min,1));
		}
		return rtrn;
	}
	
	/**
	 * Get list of all bins by specifying number of bins
	 * @param numBins Number of bins
	 * @return List of bins
	 */
	public ArrayList<ExpressionQuantileBin> getQuantileBinsByNumBins(int numBins) {
		if(numBins < 1) {
			throw new IllegalArgumentException("Number of bins must be >= 1");
		}
		ArrayList<ExpressionQuantileBin> rtrn = new ArrayList<ExpressionQuantileBin>();
		double binSize = (double) 1 / (double) numBins;
		for(int i=0; i<numBins; i++) {
			rtrn.add(new ExpressionQuantileBin(i * binSize, (i + 1) * binSize));
		}
		return rtrn;
	}
	
	private class StatisticAverageDistanceBetweenFeatures extends AbstractStatistic {

		@SuppressWarnings("synthetic-access")
		public StatisticAverageDistanceBetweenFeatures() {
			// TODO Auto-generated constructor stub
		}

		@Override
		public void calculate() throws IOException {
			makeFeatureStatsDistributions();
		}

		@Override
		public Map<ExpressionQuantileBin, EmpiricalDistribution> getDistributions() {
			return distAvgDistBetweenFeatures;
		}

		@Override
		public String statisticName() {
			return "avg_distance_between_features";
		}
		
	}
	
	private class StatisticFeatureSize extends AbstractStatistic {

		@SuppressWarnings("synthetic-access")
		public StatisticFeatureSize() {
			// TODO Auto-generated constructor stub
		}

		@Override
		public void calculate() throws IOException {
			makeDistFeatureSize();
		}

		@Override
		public Map<ExpressionQuantileBin, EmpiricalDistribution> getDistributions() {
			if(distFeatureSize == null) {
				return null;
			}
			Map<ExpressionQuantileBin, EmpiricalDistribution> rtrn = new HashMap<ExpressionQuantileBin, EmpiricalDistribution>();
			rtrn.put(allQuantiles(), distFeatureSize);
			return rtrn;
		}

		@Override
		public String statisticName() {
			return "feature_size";
		}
		
	}
	
	private class StatisticMaxDistanceBetweenConsecutiveFeatures extends AbstractStatistic {

		@SuppressWarnings("synthetic-access")
		public StatisticMaxDistanceBetweenConsecutiveFeatures() {}
		
		@Override
		public void calculate() throws IOException {
			makeConsecutiveFeatureDistributions();
		}

		@Override
		public Map<ExpressionQuantileBin, EmpiricalDistribution> getDistributions() {
			return distMaxDistanceBtwConsecutiveFeatures;
		}

		@Override
		public String statisticName() {
			return "max_distance_between_consecutive_features";
		}
		
	}
	
	private class StatisticMedianDistanceBetweenConsecutiveFeatures extends AbstractStatistic {

		@SuppressWarnings("synthetic-access")
		public StatisticMedianDistanceBetweenConsecutiveFeatures() {}
		
		@Override
		public void calculate() throws IOException {
			makeConsecutiveFeatureDistributions();
		}

		@Override
		public Map<ExpressionQuantileBin, EmpiricalDistribution> getDistributions() {
			return distMedianDistanceBtwConsecutiveFeatures;
		}

		@Override
		public String statisticName() {
			return "median_distance_between_consecutive_features";
		}
		
	}
	
	private class StatisticMedianFeatureSizeByGene extends AbstractStatistic {

		@SuppressWarnings("synthetic-access")
		public StatisticMedianFeatureSizeByGene() {}
		
		@Override
		public void calculate() throws IOException {
			makeFeatureStatsDistributions();
		}

		@Override
		public Map<ExpressionQuantileBin, EmpiricalDistribution> getDistributions() {
			return distMedianFeatureSizeByGene;
		}

		@Override
		public String statisticName() {
			return "median_feature_size_by_gene";
		}
		
	}
	
	private class StatisticNumFeaturesPerGene extends AbstractStatistic {

		@SuppressWarnings("synthetic-access")
		public StatisticNumFeaturesPerGene() {}
		
		@Override
		public void calculate() throws IOException {
			makeFeatureCoverageDistributions();
		}

		@Override
		public Map<ExpressionQuantileBin, EmpiricalDistribution> getDistributions() {
			return distNumFeaturesPerGene;
		}

		@Override
		public String statisticName() {
			return "num_features_per_gene";
		}
		
	}
	
	private class StatisticPercentFeatureCoveragePerGene extends AbstractStatistic {

		@SuppressWarnings("synthetic-access")
		public StatisticPercentFeatureCoveragePerGene() {}
		
		@Override
		public void calculate() throws IOException {
			makeFeatureCoverageDistributions();
		}

		@Override
		public Map<ExpressionQuantileBin, EmpiricalDistribution> getDistributions() {
			return distPctFeaturesPerGene;
		}

		@Override
		public String statisticName() {
			return "percent_feature_coverage_per_gene";
		}
		
	}
	
	private interface Statistic {
		
		public void calculate() throws IOException;
		public Map<ExpressionQuantileBin, EmpiricalDistribution> getDistributions();
		public String statisticName();
				
	}
	
	private abstract class AbstractStatistic implements Statistic {
		
		@Override
		public String toString() {
			return statisticName();
		}
		
	}
	
	protected ExpressionQuantileBin allQuantiles() {
		return new ExpressionQuantileBin(0,1);
	}

	
	private class ExpressionQuantileBin {
		
		private double minQuantileInclusive;
		private double maxQuantileExclusive;
		
		public ExpressionQuantileBin(double minInclusive, double maxExclusive) {
			if(minInclusive < 0) {
				throw new IllegalArgumentException("Min quantile must be >= 0");
			}
			if(maxExclusive > 1) {
				throw new IllegalArgumentException("Max quantile must be <= 1");
			}
			minQuantileInclusive = minInclusive;
			maxQuantileExclusive = maxExclusive;
		}
		
		/**
		 * Get the min of the bin (inclusive)
		 * @return The min
		 */
		public double getMinQuantileInclusive() {
			return minQuantileInclusive;
		}
		
		/**
		 * Get the max of the bin (exclusive)
		 * @return The max
		 */
		public double getMaxQuantileExclusive() {
			return maxQuantileExclusive;
		}
		
		/**
		 * Whether the quantile is in the bin
		 * @param quantile The quantile
		 * @return True if quantile is in bin or if exclusive max is 1
		 */
		public boolean inBin(double quantile) {
			if(quantile < 0 || quantile > 1) {
				throw new IllegalArgumentException("Quantile must be between 0 and 1");
			}
			if(quantile >= minQuantileInclusive && maxQuantileExclusive == 1) return true;
			return quantile >= minQuantileInclusive && quantile < maxQuantileExclusive;
		}
		
		@Override
		public boolean equals(Object o) {
			if(!o.getClass().equals(getClass())) return false;
			ExpressionQuantileBin e = (ExpressionQuantileBin)o;
			return e.getMinQuantileInclusive() == minQuantileInclusive && e.getMaxQuantileExclusive() == maxQuantileExclusive;
		}
		
		@Override
		public int hashCode() {
			return toString().hashCode();
		}
		
		@Override
		public String toString() {
			return minQuantileInclusive + "-" + maxQuantileExclusive;
		}
		
	}
	
	private class ExpressionScoreFile {
		
		private String fileName;
		
		protected ExpressionScoreFile() {
			fileName = getExpressionFileName(bedAnnotation);
		}
		
		protected boolean readFromFile() throws IOException {
			File f = new File(fileName);
			if(!f.exists()) {
				return false;
			}
			logger.info("");
			logger.info("Trying to read expression scores from file " + fileName + "...");
			FileReader r = new FileReader(fileName);
			BufferedReader b = new BufferedReader(r);
			Map<Gene, Double> levels = new TreeMap<Gene, Double>();
			Collection<String> geneNamesDone = new TreeSet<String>();
			StringParser p = new StringParser();
			while(b.ready()) {
				String line = b.readLine();
				p.parse(line);
				if(p.getFieldCount() != 2) {
					continue;
				}
				String geneName = p.asString(0);
				double score = p.asDouble(1);
				if(!genesByName.containsKey(geneName)) {
					continue;
				}
				Gene gene = genesByName.get(geneName);
				levels.put(gene, Double.valueOf(score));
				geneNamesDone.add(gene.getName());
			}
			b.close();
			if(geneNamesDone.equals(genesByName.keySet())) {
				expressionLevels = new TreeMap<Gene, Double>();
				expressionLevels.putAll(levels);
				ArrayList<Double> levelsList = new ArrayList<Double>();
				for(Gene gene : expressionLevels.keySet()) {
					levelsList.add(expressionLevels.get(gene));
				}
				sortedExpressionLevels = new double[levelsList.size()];
				for(int i=0; i < sortedExpressionLevels.length; i++) {
					sortedExpressionLevels[i] = levelsList.get(i).doubleValue();
				}
				Arrays.sort(sortedExpressionLevels);
				logger.info("Read from file.");
				return true;
			}
			logger.info("Could not read from file.");
			return false;
		}
		
		protected void writeToFile() throws IOException {
			logger.info("");
			logger.info("Writing expression scores to file " + fileName + "...");
			FileWriter w = new FileWriter(fileName);
			for(Gene gene : expressionLevels.keySet()) {
				w.write(gene.getName() + "\t" + expressionLevels.get(gene).doubleValue() + "\n");
			}
			w.close();
			logger.info("Done writing to file.");
		}
		
		private String getExpressionFileName(String bedFile) {
			return bedFile + ".expressionScores";
		}
		
	}

}
