/**
 * 
 */
package util.programs.counts;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

import org.apache.log4j.Logger;

import broad.core.overlaputils.GeneSetIntersect;
import broad.pda.annotation.BEDFileParser;

import nextgen.core.annotation.Gene;
import nextgen.core.coordinatesystem.TranscriptomeSpace;
import nextgen.core.model.ScanStatisticDataAlignmentModel;
import nextgen.core.model.score.ScanStatisticScore;

/**
 * @author prussell
 *
 */
public class ComparativeFeatureEnrichment {
	
	private Map<String, Collection<Gene>> genes;
	private Map<String, Map<String, Collection<Gene>>> features;
	private Map<String, ScanStatisticDataAlignmentModel> data;
	private TranscriptomeSpace transcriptomeSpace;
	private Map<String, Map<String, Map<Gene, Collection<Gene>>>> overlaps;
	private Map<String, Map<Gene, Collection<ScanStatisticScore>>> featureScores;
	private Map<Gene, Double> geneCounts;
	
	private static String BACKGROUND_SAMPLE_NAME = "background";
	private static Logger logger = Logger.getLogger(ComparativeFeatureEnrichment.class.getName());
	
	private ComparativeFeatureEnrichment(String backgroundBamFile, String bedGeneAnnotation, Map<String, String> bamFilesBySampleName, Map<String, String> featureBedFilesBySampleName) throws IOException {
		if(!bamFilesBySampleName.keySet().equals(featureBedFilesBySampleName.keySet())) {
			throw new IllegalArgumentException("Sample names must be the same for bam files and feature bed files.");
		}
		logger.info("");
		logger.info("Loading genes...");
		genes = BEDFileParser.loadDataByChr(new File(bedGeneAnnotation));
		transcriptomeSpace = new TranscriptomeSpace(genes);
		logger.info("");
		logger.info("Loading background alignment data...");
		geneCounts = new TreeMap<Gene, Double>();
		loadData(backgroundBamFile, bamFilesBySampleName);
		loadFeatures(featureBedFilesBySampleName);
		mapGenesToFeatures();
		scoreFeatures();
	}
	
	private void mapGenesToFeatures() {
		logger.info("");
		logger.info("Identifying each feature with overlapping genes...");
		overlaps = new TreeMap<String, Map<String, Map<Gene, Collection<Gene>>>>();
		for(String sampleName : features.keySet()) {
			logger.info(sampleName);
			overlaps.put(sampleName, GeneSetIntersect.mapGenesToOverlappers(features.get(sampleName), genes));
		}
	}
	
	private void loadData(String backgroundBamFile, Map<String, String> bamFilesBySampleName) {
		logger.info("");
		logger.info("Loading alignment data...");
		data = new TreeMap<String, ScanStatisticDataAlignmentModel>();
		for(String sampleName : bamFilesBySampleName.keySet()) {
			logger.info(sampleName);
			String bamFile = bamFilesBySampleName.get(sampleName);
			data.put(sampleName, new ScanStatisticDataAlignmentModel(bamFile, transcriptomeSpace));
		}
		logger.info(BACKGROUND_SAMPLE_NAME);
		data.put(BACKGROUND_SAMPLE_NAME, new ScanStatisticDataAlignmentModel(backgroundBamFile, transcriptomeSpace));
	}
	
	private void loadFeatures(Map<String, String> featureBedFilesBySampleName) throws IOException {
		logger.info("");
		logger.info("Loading features...");
		features = new TreeMap<String, Map<String, Collection<Gene>>>();
		for(String sampleName : featureBedFilesBySampleName.keySet()) {
			logger.info(sampleName);
			String bedFile = featureBedFilesBySampleName.get(sampleName);
			features.put(sampleName, BEDFileParser.loadDataByChr(new File(bedFile)));
		}		
	}
	
	private void scoreFeatures() {
		logger.info("");
		logger.info("Scoring features...");
		featureScores = new TreeMap<String, Map<Gene, Collection<ScanStatisticScore>>>();
		for(String sampleName : features.keySet()) {
			Map<Gene, Collection<ScanStatisticScore>> scores = new TreeMap<Gene, Collection<ScanStatisticScore>>();
			for(String chr : features.get(sampleName).keySet()) {
				logger.info(sampleName + "\t" + chr);
				for(Gene feature : features.get(sampleName).get(chr)) {
					Collection<ScanStatisticScore> scoresThisGene = new ArrayList<ScanStatisticScore>();
					for(Gene overlapper : overlaps.get(sampleName).get(chr).get(feature)) {
						double geneCount = getGeneCount(sampleName, overlapper);
						double geneLength = overlapper.getSize();
						ScanStatisticScore score = new ScanStatisticScore(data.get(sampleName), feature, geneCount, geneLength, geneCount, geneLength, false);
						scoresThisGene.add(score);
					}
					scores.put(feature, scoresThisGene);
				}
			}
			featureScores.put(sampleName, scores);
		}
	}
	
	private double getGeneCount(String sampleName, Gene gene) {
		if(geneCounts.containsKey(gene)) {
			return geneCounts.get(gene).doubleValue();
		}
		double rtrn = data.get(sampleName).getCount(gene);
		geneCounts.put(gene, Double.valueOf(rtrn));
		return rtrn;
	}
	

}
