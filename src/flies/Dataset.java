/**
 * 
 */
package flies;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.apache.commons.math3.distribution.ExponentialDistribution;

import broad.core.math.EmpiricalDistribution;
import broad.core.math.Statistics;
import broad.core.parser.CommandLineParser;
import broad.core.parser.StringParser;


/**
 * @author prussell
 * A collection of experiments
 */
public class Dataset {

	/*
	 * Column headers in data file
	 */
	private static String EXPERIMENT_NAME_COL = "Sample_ID";
	private static String EXPERIMENT_TOTAL_TIME_COL = "Total_time";
	private static String SITE_ID_COL = "Site_ID";
	private static String REL_TIME_COL = "Relative_time";
	private static String CENTROID_X_COL = "Centroid_X";
	private static String CENTROID_Y_COL = "Centroid_Y";
	
	/*
	 * Experiments by experiment name
	 */
	private Map<String, Experiment> experimentsByName;
	
	/*
	 * Output directory for analysis
	 */
	private String outDir;
	
	/*
	 * Significance cutoff
	 */
	double alpha;
	
	/*
	 * Time interval or window size
	 */
	double timeInterval;
	
	/*
	 * Max time between consecutive events in a cluster
	 */
	double clusterMaxInterEventTime;
	
	/*
	 * Step size for windows
	 */
	double stepSize;
	
	/*
	 * Number of random sites to use for permutation tests
	 */
	int numRandomPermutations;
	
	/*
	 * Instantiate with data file
	 */
	private Dataset(String dataFile, String outDirectory, double sigCutoff, double window, double step, int numPermutations, double clusterInterEventTimeCutoff) throws IOException {
		experimentsByName = new TreeMap<String, Experiment>();
		outDir = outDirectory;
		alpha = sigCutoff;
		timeInterval = window;
		stepSize = step;
		numRandomPermutations = numPermutations;
		clusterMaxInterEventTime = clusterInterEventTimeCutoff;
		File outDirFile = new File(outDir);
		@SuppressWarnings("unused")
		boolean madeDir = outDirFile.mkdir();
		parseDataFile(dataFile);
	}
	
	/*
	 * Get number of experiments
	 */
	@SuppressWarnings("unused")
	private int getNumExperiments() {
		return experimentsByName.size();
	}
	
	/*
	 * Get total number of sites in all experiments
	 */
	@SuppressWarnings("unused")
	private int getNumSites() {
		int rtrn = 0;
		for(String expName : experimentsByName.keySet()) {
			rtrn += experimentsByName.get(expName).numSites();
		}
		return rtrn;
	}
	
	/*
	 * X coordinate of centroid of the sites
	 */
	static double centroidXcoord(Collection<Site> sites) {
		double sumXcoords = 0;
		for(Site site : sites) {
			sumXcoords += site.getXCoord();
		}
		return sumXcoords / sites.size();
	}

	/*
	 * Y coordinate of centroid of the sites
	 */
	static double centroidYcoord(Collection<Site> sites) {
		double sumYcoords = 0;
		for(Site site : sites) {
			sumYcoords += site.getYCoord();
		}
		return sumYcoords / sites.size();
	}

	/*
	 * Distance from centroid to furthest site
	 */
	double radius(Collection<Site> sites) {
		double centroidX = centroidXcoord(sites);
		double centroidY = centroidYcoord(sites);
		EuclidianDistance d = new EuclidianDistance();
		double rtrn = 0;
		for(Site site : sites) {
			double dist = d.getDistance(centroidX, centroidY, site);
			if(dist > rtrn) rtrn = dist;
		}
		return rtrn;
	}

	/*
	 * For an experiment, get all sites within a radius of a point according to a distance
	 */
	static Collection<Site> getAllSitesWithinRadius(Experiment exp, double xCoord, double yCoord, double radius, Distance dist) {
		Collection<Site> rtrn = new TreeSet<Site>();
		for(Site site : exp.getSites()) {
			if(dist.getDistance(xCoord, yCoord, site) <= radius) {
				rtrn.add(site);
			}
		}			
		return rtrn;
	}

	
	/*
	 * Get total number of events in all experiments
	 */
	@SuppressWarnings("unused")
	private int getNumEvents() {
		int rtrn = 0;
		for(String expName : experimentsByName.keySet()) {
			rtrn += experimentsByName.get(expName).numEvents();
		}
		return rtrn;
	}
	
	/*
	 * Load data from file
	 */
	private void parseDataFile(String fileName) throws IOException {
		
		FileReader r = new FileReader(fileName);
		BufferedReader b = new BufferedReader(r);
		String header = b.readLine();
		StringParser p = new StringParser();
		p.parse(header);
		int expNameCol = p.getIndexFor(EXPERIMENT_NAME_COL);
		int expTotalTimeCol = p.getIndexFor(EXPERIMENT_TOTAL_TIME_COL);
		int siteIdCol = p.getIndexFor(SITE_ID_COL);
		int relTimeCol = p.getIndexFor(REL_TIME_COL);
		int xCoordCol = p.getIndexFor(CENTROID_X_COL);
		int yCoordCol = p.getIndexFor(CENTROID_Y_COL);
		
		while(b.ready()) {
			String line = b.readLine();
			p.parse(line);
			if(p.getFieldCount() == 0) continue;
			String expName = p.asString(expNameCol);
			String siteId = p.asString(siteIdCol);
			double expTotalTime = p.asDouble(expTotalTimeCol);
			double relTime = p.asDouble(relTimeCol);
			double xCoord = p.asDouble(xCoordCol);
			double yCoord = p.asDouble(yCoordCol);
			if(!experimentsByName.containsKey(expName)) {
				Experiment exp = new Experiment(expName, expTotalTime);
				Site site = new Site(exp, siteId, xCoord, yCoord);
				site.addEvent(relTime);
				exp.addSite(site);
				experimentsByName.put(expName, exp);
				continue;
			}
			Experiment exp = experimentsByName.get(expName);
			boolean hasSite = false;
			for(Site site : exp.getSites()) {
				if(site.getSiteId().equals(siteId)) {
					site.addEvent(relTime);
					hasSite = true;
					continue;
				}
			}
			if(!hasSite) {
				Site site = new Site(exp, siteId, xCoord, yCoord);
				site.addEvent(relTime);
				exp.addSite(site);
				
			}
		}
		
		r.close();
		b.close();
		
	}
	
	/*
	 * Distribution of number of events following an event at the same site within the time interval
	 * All sites counted separately
	 * All experiments
	 */
	private double[] countDistributionConsecutiveEvents() {
		// Get the max possible count
		int maxCount = 0;
		for(String name : experimentsByName.keySet()) {
			Experiment exp = experimentsByName.get(name);
			for(Site site : exp.getSites()) {
				for(Event event : site.getEvents()) {
					int count = site.countEventsAfter(event);
					if(count > maxCount) {
						maxCount = count;
					}
				}
			}
		}
		// Populate the counts distribution
		int[] counts = new int[maxCount + 1];
		for(int i=0; i < counts.length; i++) counts[i] = 0;
		for(String name : experimentsByName.keySet()) {
			Experiment exp = experimentsByName.get(name);
			for(Site site : exp.getSites()) {
				for(Event event : site.getEvents()) {
					int count = site.countEventsAfter(event);
					counts[count]++;
				}
			}
		}
		// Normalize to sum of 1
		int total = 0;
		for(int i=0; i < counts.length; i++) total += counts[i];
		double[] rtrn = new double[counts.length];
		for(int i=0; i < rtrn.length; i++) rtrn[i] = (double)counts[i]/(double)total;
		return rtrn;
	}

	
	/*
	 * For randomized permutations of all experiments, get distribution of number of events following an event within the time interval
	 */
	static double[] countDistributionConsecutiveEventsRandomized(ArrayList<Experiment> randomPermutations) {
		// Get the max possible count
		int maxCount = 0;
		for(Experiment randExp : randomPermutations) {
			for(Site site : randExp.getSites()) {
				for(Event event : site.getEvents()) {
					int count = site.countEventsAfter(event);
					if(count > maxCount) {
						maxCount = count;
					}
				}
			}
		}
		// Populate the counts distribution
		int[] counts = new int[maxCount + 1];
		for(int i=0; i < counts.length; i++) counts[i] = 0;
		for(Experiment randExp : randomPermutations) {
			for(Site site : randExp.getSites()) {
				for(Event event : site.getEvents()) {
					int count = site.countEventsAfter(event);
					counts[count]++;
				}
			}
		}
		// Normalize to sum of 1
		int total = 0;
		for(int i=0; i < counts.length; i++) total += counts[i];
		double[] rtrn = new double[counts.length];
		for(int i=0; i < rtrn.length; i++) rtrn[i] = (double)counts[i]/(double)total;
		return rtrn;
		
	}

	/**
	 * Across all the experiments in the set, the percentage of events that are followed by another event at the same site within the time span
	 * @param experiments The experiments
	 * @param timeSpan The time span
	 * @return The percentage of events followed by another event at the site within the time span
	 */
	private static double consecutiveEventPercentage(Collection<Experiment> experiments, double timeSpan) {
		double numWithNext = 0;
		double totalEvents = 0;
		for(Experiment exp : experiments) {
			for(Site site : exp.getSites()) {
				for(Event event : site.getEvents()) {
					totalEvents++;
					try {
						double time = site.timeToNextEvent(event);
						if(time <= timeSpan) numWithNext++;
					} catch(IllegalArgumentException e) {
						continue;
					}
				}
			}
		}
		return numWithNext / totalEvents;
	}

	
	/**
	 * The percentage of events that are followed by another event at the same site within the time span
	 * @param exp The experiment
	 * @param timeSpan The time span
	 * @return The percentage of events followed by another event at the site within the time span
	 */
	private static double consecutiveEventPercentage(Experiment exp, double timeSpan) {
		ArrayList<Experiment> experiments = new ArrayList<Experiment>();
		experiments.add(exp);
		return consecutiveEventPercentage(experiments, timeSpan);
	}

	
	/**
	 * Write a distribution to a file
	 * @param outFile Output file
	 * @param distribution The distribution where element i is number of data points in bin i
	 * @throws IOException
	 */
	private static void writeDistribution(String outFile, double[] distribution) throws IOException {
		FileWriter w = new FileWriter(outFile);
		for(int i=0; i < distribution.length; i++) {
			w.write(i + "\t" + distribution[i] + "\n");
		}
		w.close();
	}
	
	/*
	 * Write list of all event clusters to file
	 */
	private void writeEventClusterList() throws IOException {
		
		String outFile = outDir + "/event_cluster_list_" + clusterMaxInterEventTime;
		
		FileWriter w = new FileWriter(outFile);
		for(String expName : experimentsByName.keySet()) {
			Experiment exp = experimentsByName.get(expName);
			Map<Site, Collection<Site.EventCluster>> clusters = exp.findEventClusters();
			for(Site site : clusters.keySet()) {
				for(Site.EventCluster cluster : clusters.get(site)) {
					w.write(cluster.toString() + "\n");
				}
			}
		}
		w.close();
	}
	
	/*
	 * Write the percentage of events having another event following within the time interval
	 */
	private void writeConsecutiveEventsPercentage() throws IOException {
		String outConsecutiveFile = outDir + "/consecutive_events_percentage_" + timeInterval;
		System.out.println("Writing percentage of events with a consecutive event to file " + outConsecutiveFile);
		FileWriter w = new FileWriter(outConsecutiveFile);
		String header = "Experiment_name\tPct_events_with_consecutive\tRandomized_times_pct_events_with_consecutive\tPct_increase_real_over_random\n";
		w.write(header);
		Collection<Experiment> allExperiments = new ArrayList<Experiment>();
		Collection<Experiment> allExperimentsRandomized = new ArrayList<Experiment>();
		for(String expName : experimentsByName.keySet()) {
			String line = expName + "\t";
			Experiment exp = experimentsByName.get(expName);
			allExperiments.add(exp);
			double realPct = consecutiveEventPercentage(exp, timeInterval);
			ArrayList<Experiment> rand = exp.getRandomPermutations();
			allExperimentsRandomized.addAll(rand);
			double randomizedPct = consecutiveEventPercentage(rand, timeInterval);
			double pctIncrease = (realPct - randomizedPct) / randomizedPct;
			line += realPct + "\t" + randomizedPct + "\t" + pctIncrease + "\n";
			w.write(line);
		}
		String line = "all_experiments\t";
		double realPct = consecutiveEventPercentage(allExperiments, timeInterval);
		double randomizedPct = consecutiveEventPercentage(allExperimentsRandomized, timeInterval);
		double pctIncrease = (realPct - randomizedPct) / randomizedPct;
		line += realPct + "\t" + randomizedPct + "\t" + pctIncrease + "\n";
		w.write(line);
		w.close();
	}
	
	/*
	 * Write distribution of number of events following another event
	 */
	private void writeConsecutiveEventsDistribution() throws IOException {
		String outConsecutiveFileAll = outDir + "/consecutive_events_distribution_all_experiments_" + timeInterval;
		System.out.println("Writing consecutive events distribution to file " + outConsecutiveFileAll);
		writeDistribution(outConsecutiveFileAll, countDistributionConsecutiveEvents());		
		for(String expName : experimentsByName.keySet()) {
			Experiment exp = experimentsByName.get(expName);
			String outConsecutiveFile = outDir + "/consecutive_events_distribution_" + expName + "_" + timeInterval;
			System.out.println("Writing consecutive events distribution to file " + outConsecutiveFile);
			writeDistribution(outConsecutiveFile, exp.countDistributionConsecutiveEvents());
		}
	}

	/*
	 * Write distribution of number of events following another event among random permutations of the experiment
	 */
	private void writeConsecutiveEventsDistributionRandomized() throws IOException {
		String outConsecutiveFileAll = outDir + "/consecutive_events_distribution_randomized_experiments_all_" + timeInterval;
		System.out.println("Writing consecutive events distribution for randomized experiments to file " + outConsecutiveFileAll);
		System.out.println("Writing consecutive events distribution for randomized experiments to file " + outConsecutiveFileAll);
		ArrayList<Experiment> randExperiments = new ArrayList<Experiment>();
		for(String expName : experimentsByName.keySet()) {
			randExperiments.addAll(experimentsByName.get(expName).getRandomPermutations());
		}
		writeDistribution(outConsecutiveFileAll, countDistributionConsecutiveEventsRandomized(randExperiments));
		for(String expName : experimentsByName.keySet()) {
			Experiment exp = experimentsByName.get(expName);
			String outConsecutiveFile = outDir + "/consecutive_events_distribution_randomized_experiments_" + expName + "_" + timeInterval;
			System.out.println("Writing consecutive events distribution for randomized experiments to file " + outConsecutiveFile);
			writeDistribution(outConsecutiveFile, exp.countDistributionConsecutiveEventsRandomized());
		}
	}

	/*
	 * Distribution of number of events per window
	 */
	private void writeEventsPerWindowDistribution() throws IOException {
		for(String expName : experimentsByName.keySet()) {
			Experiment exp = experimentsByName.get(expName);
			String outWindowFile = outDir + "/events_per_window_" + expName + "_" + timeInterval + "_" + stepSize;
			System.out.println("Writing distribution of number of events per window to file " + outWindowFile);
			writeDistribution(outWindowFile, exp.countDistributionEventsPerWindow());
		}
	}

	/*
	 * Write site wide significant clusters
	 */
	private void writeSiteWideSigClusters() throws IOException {
		String outSiteWideSigClusterFile = outDir + "/sig_clusters_site_wide_" + clusterMaxInterEventTime + "_" + alpha;
		FileWriter writeSiteWideSigClusters = new FileWriter(outSiteWideSigClusterFile);
		System.out.println("Writing site-wide significant clusters to file " + outSiteWideSigClusterFile);
		for(String expName : experimentsByName.keySet()) {
			Experiment exp = experimentsByName.get(expName);
			for(Site site : exp.getSites()) {
				// Write Site wide significant clusters
				ArrayList<Site.EventCluster> siteWideSigClusters = site.significantEventClusters();
				for(Site.EventCluster cluster : siteWideSigClusters)	{
					writeSiteWideSigClusters.write(cluster.toString() + "\n");
				}
			}			
		}
		writeSiteWideSigClusters.close();
	}
	
	/*
	 * Write tables of site scores to files for individual experiments and all experiments combined
	 */
	private void writeSiteScores() throws IOException {
		String outSiteScoreFile = outDir + "/site_scores_all_experiments";
		FileWriter siteScoreWriter = new FileWriter(outSiteScoreFile);
		System.out.println("Writing site scores (experiments combined) to file " + outSiteScoreFile);
		String header = "Site_ID\t";
		header += "Avg_time_per_event\t";
		header += "Cluster_pval\t";
		header += "Cluster_score\t";
		siteScoreWriter.write(header + "\n");
		for(String expName : experimentsByName.keySet()) {
			String outScoreFileExp = outDir + "/site_scores_" + expName;
			FileWriter siteScoreWriterExp = new FileWriter(outScoreFileExp);
			System.out.println("Writing site scores for experiment " + expName + " to file " + outScoreFileExp);
			siteScoreWriterExp.write(header + "\n");
			Experiment exp = experimentsByName.get(expName);
			Map<Site, SiteScore> scores = exp.siteScores();
			for(Site site : scores.keySet()) {
				String line = site.getSiteId() + "\t";
				SiteScore score = scores.get(site);
				line += score.getAvgTimePerEvent() + "\t";
				line += score.getSiteClusterPval() + "\t";
				line += score.getSiteClusterScore() + "\t";
				siteScoreWriter.write(line + "\n");
				siteScoreWriterExp.write(line + "\n");
			}
			siteScoreWriterExp.close();
		}
		siteScoreWriter.close();
	}
	
	/*
	 * Write experiment wide significant clusters
	 */
	private void writeExpWideSigClusters() throws IOException {
		String outExpWideSigClusterFile = outDir + "/sig_clusters_exp_wide_" + clusterMaxInterEventTime + "_" + alpha;
		FileWriter writeExpWideSigClusters = new FileWriter(outExpWideSigClusterFile);
		System.out.println("Writing experiment-wide significant clusters to file " + outExpWideSigClusterFile);
		for(String expName : experimentsByName.keySet()) {
			Experiment exp = experimentsByName.get(expName);
			Map<Integer, Double> expWideScoreCutoffs = exp.scoreCutoffs();
			for(Site site : exp.getSites()) {
				// Write experiment wide significant clusters
				ArrayList<Site.EventCluster> expWideSigClusters = site.significantEventClusters(expWideScoreCutoffs);
				for(Site.EventCluster cluster : expWideSigClusters)	{
					writeExpWideSigClusters.write(cluster.toString() + "\n");
				}

			}
		}
		writeExpWideSigClusters.close();
	}

	/*
	 * Write summary of neighborhoods to file
	 */
	@SuppressWarnings("unused")
	private void writeNeighborhoodSummary() throws IOException {
		String outNeighborhoodSummaryFile = outDir + "/neighborhood_summary";
		FileWriter w = new FileWriter(outNeighborhoodSummaryFile);
		System.out.println("Writing neighborhood summary to file " + outNeighborhoodSummaryFile);
		for(String exp : experimentsByName.keySet()) {
			for(Neighborhood n : experimentsByName.get(exp).getAllNeighborhoods()) {
				w.write(n.toString() + "\n");
			}
		}
		w.close();
	}

	/*
	 * Write summary of event clusters to file
	 */
	private void writeEventClusterSummary() throws IOException {
		
		String outFile = outDir + "/event_cluster_summary_" + clusterMaxInterEventTime;
		
		String longestClusterHeader = "Experiment_name\t";
		longestClusterHeader += "Total_sites\t";
		longestClusterHeader += "Num_sites_no_clusters\t";
		longestClusterHeader += "Num_sites_longest_cluster_2\t";
		longestClusterHeader += "Num_sites_longest_cluster_3\t";
		longestClusterHeader += "Num_sites_longest_cluster_4\t";
		longestClusterHeader += "Num_sites_longest_cluster_5\t";
		longestClusterHeader += "Num_sites_longest_cluster_>5\n";
		
		String numClustersHeader = "Experiment_name\t";
		numClustersHeader += "Total_sites\t";
		numClustersHeader += "Num_sites_no_clusters\t";
		numClustersHeader += "Num_sites_1_cluster\t";
		numClustersHeader += "Num_sites_2_clusters\t";
		numClustersHeader += "Num_sites_3_clusters\t";
		numClustersHeader += "Num_sites_>3_clusters\n";
		

		Map<String, String> linesToWriteLongestCluster = new TreeMap<String, String>();
		Map<String, String> linesToWriteNumClusters = new TreeMap<String, String>();
		
		for(String expName : experimentsByName.keySet()) {
			
			Experiment exp = experimentsByName.get(expName);
			TreeMap<Integer,Integer> countBySize = new TreeMap<Integer,Integer>();
			TreeMap<Integer,Integer> countByNumClusters = new TreeMap<Integer, Integer>();
			
			for(Site site : exp.getSites()) {
				
				Integer largestCluster = Integer.valueOf(site.largestClusterSizeByInterEventTime());
				if(!countBySize.containsKey(largestCluster)) {
					countBySize.put(largestCluster, Integer.valueOf(1));
				}
				else {
					int oldCount = countBySize.get(largestCluster).intValue();
					countBySize.put(largestCluster, Integer.valueOf(oldCount + 1));
				}
				
				
				Integer numClusters = Integer.valueOf(site.findEventClustersByInterEventTime().size());
				if(!countByNumClusters.containsKey(numClusters)) {
					countByNumClusters.put(numClusters, Integer.valueOf(1));
				}
				else {
					int oldCount2 = countByNumClusters.get(numClusters).intValue();
					countByNumClusters.put(numClusters, Integer.valueOf(oldCount2 + 1));
				}
			}
			
			String longestClusterLine = exp.getName() + "\t";
			longestClusterLine += exp.numSites() + "\t";
			longestClusterLine += (countBySize.containsKey(Integer.valueOf(0)) ? countBySize.get(Integer.valueOf(0)) : "0") + "\t";
			longestClusterLine += (countBySize.containsKey(Integer.valueOf(2)) ? countBySize.get(Integer.valueOf(2)) : "0") + "\t";
			longestClusterLine += (countBySize.containsKey(Integer.valueOf(3)) ? countBySize.get(Integer.valueOf(3)) : "0") + "\t";
			longestClusterLine += (countBySize.containsKey(Integer.valueOf(4)) ? countBySize.get(Integer.valueOf(4)) : "0") + "\t";
			longestClusterLine += (countBySize.containsKey(Integer.valueOf(5)) ? countBySize.get(Integer.valueOf(5)) : "0") + "\t";
			Map<Integer, Integer> remainder = countBySize.tailMap(Integer.valueOf(6));
			int remainderSize = 0;
			for(Integer i : remainder.keySet()) {
				remainderSize += remainder.get(i).intValue();
			}
			longestClusterLine += remainderSize + "\n";
			
			String numClustersLine = exp.getName() + "\t";
			numClustersLine += exp.numSites() + "\t";
			numClustersLine += (countByNumClusters.containsKey(Integer.valueOf(0)) ? countByNumClusters.get(Integer.valueOf(0)) : "0") + "\t";
			numClustersLine += (countByNumClusters.containsKey(Integer.valueOf(0)) ? countByNumClusters.get(Integer.valueOf(1)) : "0") + "\t";
			numClustersLine += (countByNumClusters.containsKey(Integer.valueOf(2)) ? countByNumClusters.get(Integer.valueOf(2)) : "0") + "\t";
			numClustersLine += (countByNumClusters.containsKey(Integer.valueOf(3)) ? countByNumClusters.get(Integer.valueOf(3)) : "0") + "\t";
			int remainderSize2 = 0;
			Map<Integer, Integer> remainder2 = countByNumClusters.tailMap(Integer.valueOf(4));
			for(Integer i : remainder2.keySet()) {
				remainderSize2 += remainder2.get(i).intValue();
			}
			numClustersLine += remainderSize2 + "\n";
			
			linesToWriteLongestCluster.put(expName, longestClusterLine);
			linesToWriteNumClusters.put(expName, numClustersLine);
			
		}
		
		FileWriter w = new FileWriter(outFile);
		
		w.write(longestClusterHeader);
		for(String expName : experimentsByName.keySet()) {
			w.write(linesToWriteLongestCluster.get(expName));
		}
		w.write("\n\n");
		
		w.write(numClustersHeader);
		for(String expName : experimentsByName.keySet()) {
			w.write(linesToWriteNumClusters.get(expName));
		}

		w.close();
	}

	/*
	 * Remove any event cluster that is fully nested inside another cluster in the set
	 */
	static ArrayList<Site.EventCluster> collapseFullyNestedClusters(ArrayList<Site.EventCluster> clusters) {
		
		if(clusters.size() <= 1) return clusters;

		boolean[] keep = new boolean[clusters.size()];
		for(int i=0; i < clusters.size(); i++) keep[i] = true;
		
		for(int i=0; i < clusters.size(); i++) {
			for(int j=0; j < clusters.size(); j++) {
				if(i == j) continue;
				if(clusters.get(i).isNestedInside(clusters.get(j))) {
					keep[i] = false;
				}
			}
		}
		
		ArrayList<Site.EventCluster> rtrn = new ArrayList<Site.EventCluster>();
		for(int i=0; i < clusters.size(); i++) {
			if(keep[i]) {
				rtrn.add(clusters.get(i));
			}
		}
		
		return rtrn;
		
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {

		CommandLineParser p = new CommandLineParser();
		p.addStringArg("--datafile", "Data file", true);
		p.addStringArg("--outdir", "Output directory", true);
		p.addDoubleArg("--window", "Time window size", false, 1);
		p.addDoubleArg("--step", "Step size for window counts (seconds)", false, 1);
		p.addDoubleArg("--inter_event_cutoff", "Max inter event time for cluster", false, 30);
		p.addDoubleArg("--alpha", "Alpha for significance", false, 0.05);
		p.addIntArg("--numrand", "Number of random sites/experiments to generate for signficance", false, 10000);
		p.addBooleanArg("--exp_wide_clusters", "Find and write experiment-wide significant clusters", false, false);
		p.addBooleanArg("--site_wide_clusters", "Find and write site-wide significant clusters", false, false);
		p.addBooleanArg("--consecutive_events_distribution", "Write distribution of number of events following event", false, false);
		p.addBooleanArg("--consecutive_events_distribution_randomized", "Write distribution of number of events following event in randomized experiments", false, false);
		p.addBooleanArg("--event_cluster_list", "Write list of all event clusters", false, false);
		p.addBooleanArg("--event_cluster_summary", "Write summary of event clusters", false, false);
		p.addBooleanArg("--events_per_window_distribution", "Write distribution of number of events per window", false, false);
		p.addBooleanArg("--consecutive_events_percentage", "Write percentage of events followed by another event within time span", false, false);
		p.addBooleanArg("--site_scores", "Write table of scores of each site", false, false);
		p.parse(args);
		String fileName = p.getStringArg("--datafile");
		String outDir = p.getStringArg("--outdir");
		double timeInterval = p.getDoubleArg("--window");
		double stepSize = p.getDoubleArg("--step");
		double alpha = p.getDoubleArg("--alpha");
		int numRand = p.getIntArg("--numrand");
		double interEvent = p.getDoubleArg("--inter_event_cutoff");
		boolean expWideClusters = p.getBooleanArg("--exp_wide_clusters");
		boolean siteWideClusters = p.getBooleanArg("--site_wide_clusters");
		boolean consecutiveEventsDist = p.getBooleanArg("--consecutive_events_distribution");
		boolean consecutiveEventsDistRandomized = p.getBooleanArg("--consecutive_events_distribution_randomized");
		boolean clusterList = p.getBooleanArg("--event_cluster_list");
		boolean clusterSummary = p.getBooleanArg("--event_cluster_summary");
		boolean eventsPerWindowDist = p.getBooleanArg("--events_per_window_distribution");
		boolean consecutiveEventsPct = p.getBooleanArg("--consecutive_events_percentage");
		boolean siteScores = p.getBooleanArg("--site_scores");
		
		Dataset fd = new Dataset(fileName, outDir, alpha, timeInterval, stepSize, numRand, interEvent);

		if(consecutiveEventsDist) fd.writeConsecutiveEventsDistribution();
		if(consecutiveEventsDistRandomized) fd.writeConsecutiveEventsDistributionRandomized();
		if(clusterList) fd.writeEventClusterList();
		if(clusterSummary) fd.writeEventClusterSummary();
		if(eventsPerWindowDist) fd.writeEventsPerWindowDistribution();
		if(siteWideClusters) fd.writeSiteWideSigClusters();
		if(expWideClusters) fd.writeExpWideSigClusters();
		if(consecutiveEventsPct) fd.writeConsecutiveEventsPercentage();
		if(siteScores) fd.writeSiteScores();
		
	}
	
	/*
	 * An experiment with sites and events
	 */
	private class Experiment implements Comparable<Experiment> {
		
		/*
		 * Experiment name
		 */
		private String expName;
		
		/*
		 * Total time
		 */
		private double expTotalTime;
		
		/*
		 * Cutoff for experiment wide cluster score significance by cluster size
		 */
		private Map<Integer, Double> scoreCutoffs;
		
		/*
		 * The sites
		 */
		private Collection<Site> sites;
		
		/*
		 * Whether experiment wide score cutoffs have been calculated
		 */
		private boolean hasScoreCutoffs;
		
		/*
		 * Random permutations of the event times within each site
		 */
		private ArrayList<Experiment> randomPermutations;
		
		/*
		 * Whether the random permutations have been generated
		 */
		private boolean hasRandomPermutations;
		
		/*
		 * Construct with name and total time
		 */
		public Experiment(String name, double totalTime) {
			expName = name;
			expTotalTime = totalTime;
			sites = new TreeSet<Site>();
			hasScoreCutoffs = false;
			hasRandomPermutations = false;
		}
		
		/*
		 * (non-Javadoc)
		 * @see java.lang.Object#toString()
		 * String with helpful information
		 */
		@Override
		public String toString() {
			return "experiment_name:" + expName + ";total_time:" + Double.valueOf(expTotalTime).toString() + ";num_sites:" + Integer.valueOf(numSites()).toString();
		}
		
		/*
		 * Get experiment name
		 */
		public String getName() {
			return expName;
		}

		/*
		 * Get number of sites
		 */
		public int numSites() {
			return sites.size();
		}
		
		/*
		 * Get total number of events at all sites
		 */
		public int numEvents() {
			int rtrn = 0;
			for(Site site : sites) {
				rtrn += site.numEvents();
			}
			return rtrn;
		}
		
		/*
		 * Add a site
		 */
		public void addSite(Site site) {
			if(!sites.contains(site)) {
				sites.add(site);
			}
		}
		
		/*
		 * Get the sites
		 */
		public Collection<Site> getSites() {
			return sites;
		}
		
		/*
		 * Get scores for each site
		 */
		public Map<Site, SiteScore> siteScores() {
			Map<Site, SiteScore> rtrn = new TreeMap<Site, SiteScore>();
			for(Site site : sites) {
				SiteScore score = new SiteScore(site);
				rtrn.put(site, score);
			}
			return rtrn;
		}
		
		/*
		 * Get total time
		 */
		public double getTotalTime() {
			return expTotalTime;
		}
		
		/*
		 * Get the set of permuted experiments with event times randomized within each site
		 */
		public ArrayList<Experiment> getRandomPermutations() {
			if(!hasRandomPermutations) {
				generateRandomPermutations();
			}
			return randomPermutations;
		}
		
		/*
		 * Generate random permutations of event times within each site
		 */
		private void generateRandomPermutations() {
			System.out.println("Generating " + numRandomPermutations + " random permutations of experiment " + expName);
			randomPermutations = new ArrayList<Experiment>();
			for(int i=0; i < numRandomPermutations; i++) {
				randomPermutations.add(randomizeEventTimes());
			}
			hasRandomPermutations = true;
		}
		
		/*
		 * For each cluster size, get the 1 - alpha quantile of max cluster score for the whole experiment
		 * where event times at each site have been randomized
		 * 
		 * numRandSites random experiments are generated
		 * For each cluster size the max score of any cluster of that size over all sites is calculated
		 */
		private Map<Integer, Double> clusterScoreCutoffEntireExperiment() {
			if(!hasRandomPermutations) {
				generateRandomPermutations();
			}
			Map<Integer, Double> rtrn = new TreeMap<Integer, Double>();
			// Generate random experiments
			// Get max score for each cluster size
			for(int size = 2; size <= getMaxNumEventsPerSite(); size++) {
				double[] randScores = new double[numRandomPermutations];
				for(int i=0; i < randomPermutations.size(); i++) {
					Experiment rand = randomPermutations.get(i);
					double randScore = rand.maxClusterScore(size);
					randScores[i] = randScore;
				}			
				Arrays.sort(randScores);
				rtrn.put(Integer.valueOf(size), Double.valueOf(Statistics.quantile(randScores, 1 - alpha)));
			}
			return rtrn;
		}

		/*
		 * Get the maximum number of events at any site
		 */
		private int getMaxNumEventsPerSite() {
			int rtrn = 0;
			for(Site site : sites) {
				int numEvents = site.numEvents();
				if(numEvents > rtrn) rtrn = numEvents;
			}
			return rtrn;
		}

		/*
		 * Initialize the map of experiment wide significant cluster score cutoffs by cluster size
		 */
		private void calculateScoreCutoffs() {
			System.out.println("Calculating experiment wide significant cluster score cutoffs for experiment " + expName + "...");
			scoreCutoffs = clusterScoreCutoffEntireExperiment();
			hasScoreCutoffs = true;
		}
		
		/*
		 * Get cluster score cutoff for experiment wide significance by cluster size
		 */
		public Map<Integer, Double> scoreCutoffs() {
			if(hasScoreCutoffs) return scoreCutoffs;
			calculateScoreCutoffs();
			return scoreCutoffs;
		}
		
		/*
		 * Get a new experiment where each site has had its event times randomly shuffled within the experiment time span
		 */
		private Experiment randomizeEventTimes() {
			Experiment rtrn = new Experiment(expName, expTotalTime);
			for(Site site : sites) {
				rtrn.addSite(site.randomizeEventTimes());
			}
			return rtrn;
		}
		
		
		/*
		 * The maximal score of all possible clusters (consecutive events) of the specified size from any site
		 * Calls the score method in EventCluster
		 */
		public double maxClusterScore(int clusterSize) {
			Collection<Site.EventCluster> allClusters = getAllEventClusters(clusterSize);
			double rtrn = 0;
			for(Site.EventCluster cluster : allClusters) {
				double score = cluster.score();
				if(score > rtrn) rtrn = score;
			}
			return rtrn;
		}

		/*
		 * Get all event clusters of the specified size from any site
		 */
		private Collection<Site.EventCluster> getAllEventClusters(int clusterSize) {
			Collection<Site.EventCluster> rtrn = new ArrayList<Site.EventCluster>();
			for(Site site : sites) {
				rtrn.addAll(site.getAllEventClusters(clusterSize));
			}
			return rtrn;
		}

		/*
		 * Get event clusters for all sites
		 */
		Map<Site, Collection<Site.EventCluster>> findEventClusters() {
			Map<Site, Collection<Site.EventCluster>> rtrn = new TreeMap<Site, Collection<Site.EventCluster>>();
			for(Site site : sites) {
				Collection<Site.EventCluster> clusters = site.findEventClustersByInterEventTime();
				if(clusters.isEmpty()) continue;
				rtrn.put(site, clusters);
			}
			return rtrn;
		}
			
		
		/*
		 * Distribution of number of events following an event at the same site within the time interval
		 * All sites counted separately
		 */
		public double[] countDistributionConsecutiveEvents() {
			// Get the max possible count
			int maxCount = 0;
			for(Site site : getSites()) {
				for(Event event : site.getEvents()) {
					int count = site.countEventsAfter(event);
					if(count > maxCount) {
						maxCount = count;
					}
				}
			}
			// Populate the counts distribution
			int[] counts = new int[maxCount + 1];
			for(int i=0; i < counts.length; i++) counts[i] = 0;
			for(Site site : getSites()) {
				for(Event event : site.getEvents()) {
					int count = site.countEventsAfter(event);
					counts[count]++;
				}
			}
			// Normalize to sum of 1
			int total = 0;
			for(int i=0; i < counts.length; i++) total += counts[i];
			double[] rtrn = new double[counts.length];
			for(int i=0; i < rtrn.length; i++) rtrn[i] = (double)counts[i]/(double)total;
			return rtrn;
		}
		
		/*
		 * For randomized permutations of the experiment, get distribution of number of events following an event within the time interval
		 */
		public double[] countDistributionConsecutiveEventsRandomized() {
			if(!hasRandomPermutations) {
				generateRandomPermutations();
			}
			return Dataset.countDistributionConsecutiveEventsRandomized(randomPermutations);
		}
		
		/*
		 * Distribution of number of events in time interval
		 * All sites counted separately
		 */
		public double[] countDistributionEventsPerWindow() {
			// Get the max possible count
			double currStart = 0;
			double currEnd = 0;
			int maxCount = 0;
			while(currStart <= getTotalTime()) {
				currEnd = currStart + timeInterval;
				for(Site site : getSites()) {
					int count = site.countEventsInWindow(currStart, currEnd);
					if(count > maxCount) {
						maxCount = count;
					}
				}
				currStart += stepSize;
			}
			// Populate the counts distribution
			int[] counts = new int[maxCount + 1];
			for(int i=0; i < counts.length; i++) counts[i] = 0;
			currStart = 0;
			currEnd = 0;
			while(currStart <= getTotalTime()) {
				currEnd = currStart + timeInterval;
				for(Site site : getSites()) {
					int count = site.countEventsInWindow(currStart, currEnd);
					counts[count]++;
				}
				currStart += stepSize;
			}
			// Normalize to sum of 1
			int total = 0;
			for(int i=0; i < counts.length; i++) total += counts[i];
			double[] rtrn = new double[counts.length];
			for(int i=0; i < rtrn.length; i++) rtrn[i] = (double)counts[i]/(double)total;
			return rtrn;			
		}
		
		/*
		 * Get all neighborhoods (subsets of sites contained within a radius of a center point)
		 */
		public Collection<Neighborhood> getAllNeighborhoods() {
			throw new UnsupportedOperationException("TODO");
			
			//TODO
		}
		
		/*
		 * (non-Javadoc)
		 * @see java.lang.Object#equals(java.lang.Object)
		 * Equal iff experiment name, total time and set of sites are equal
		 */
		@Override
		public boolean equals(Object o) {
			Experiment e = (Experiment)o;
			if(!getName().equals(e.getName())) return false;
			if(getTotalTime() != e.getTotalTime()) return false;
			return getSites().equals(e.getSites());
			
		}

		/*
		 * (non-Javadoc)
		 * @see java.lang.Comparable#compareTo(java.lang.Object)
		 * Compare experiment name, then total time, then number of sites
		 */
		@Override
		public int compareTo(Experiment other) {
			if(!getName().equals(other.getName())) {
				return getName().compareTo(other.getName());
			}
			if(getTotalTime() != other.getTotalTime()) {
				return Double.valueOf(getTotalTime()).compareTo(Double.valueOf(other.getTotalTime()));
			}
			if(numSites() != other.numSites()) {
				return Double.valueOf(numSites()).compareTo(Double.valueOf(other.numSites()));
			}
			throw new IllegalArgumentException("Can't compare experiments with same name, same total time and same number of sites.");
		}

		/*
		 * (non-Javadoc)
		 * @see java.lang.Object#hashCode()
		 * The hash code of toString() representation
		 */
		@Override
		public int hashCode() {
			return toString().hashCode();
		}

		
	}
	
	/*
	 * A site
	 */
	private class Site implements Comparable<Site>{
		
		/*
		 * The experiment that the site belongs to
		 */
		private Experiment siteExp;
		
		/*
		 * Whether site wide score cutoffs have been calculated
		 */
		private boolean hasScoreCutoffs;

		/*
		 * The site ID
		 */
		String siteID;
		
		/*
		 * The collection of events at the site
		 */
		private TreeSet<Event> siteEvents;
		
		/*
		 * X coordinate of centroid
		 */
		private double xCoord;
		
		/*
		 * Y coordinate of centroid
		 */
		private double yCoord;
		
		/*
		 * Cutoff for site wide cluster score significance by cluster size
		 */
		private Map<Integer, Double> scoreCutoffs;

		
		/*
		 * Construct with parent experiment and site ID
		 */
		public Site(Experiment experiment, String siteId, double centroidXCoord, double centroidYCoord) {
			siteExp = experiment;
			siteID = siteId;
			xCoord = centroidXCoord;
			yCoord = centroidYCoord;
			siteEvents = new TreeSet<Event>();
			hasScoreCutoffs = false;
		}
				
		/*
		 * Get a copy of the site
		 */
		@SuppressWarnings("unused")
		public Site copy() {
			Site rtrn = new Site(siteExp, siteID, xCoord, yCoord);
			for(Event event : siteEvents) {
				rtrn.addEvent(event.getRelTime());
			}
			return rtrn;
		}
		
		/*
		 * (non-Javadoc)
		 * @see java.lang.Object#toString()
		 * String with helpful information
		 */
		@Override
		public String toString() {
			return "experiment_name:" + siteExp.getName() + ";site_name:" + siteID + ";num_events:" + Integer.valueOf(numEvents()).toString();
		}
		
		/*
		 * Get the parent experiment
		 */
		public Experiment getExperiment() {
			return siteExp;
		}
		
		/*
		 * Get the site ID
		 */
		public String getSiteId() {
			return siteID;
		}
		
		/*
		 * Get the number of events at the site
		 */
		public int numEvents() {
			return siteEvents.size();
		}
		
		/*
		 * Add an event
		 */
		public void addEvent(double relTime) {
			siteEvents.add(new Event(siteExp, this, relTime));
		}
		
		/*
		 * Get all the event times
		 */
		public TreeSet<Double> getEventTimes() {
			return getEventTimes(-1, getExperiment().getTotalTime() + 1, true);
		}
		
		/*
		 * Get set of times between consecutive events
		 * The set is sorted by inter-event time
		 * Includes an entry for time to first event + time after last event
		 * If there is only one event, contains a single entry which is the total experiment time
		 * Throws exception if there are no events
		 */
		public TreeSet<Double> interEventTimes() {
			
			if(numEvents() == 0) {
				throw new IllegalStateException("Site has no events");
			}
			
			TreeSet<Double> eventTimes = getEventTimes();
			TreeSet<Double> rtrn = new TreeSet<Double>();
			Iterator<Double> iter = eventTimes.iterator();
			double curr = iter.next().doubleValue();
			
			while(iter.hasNext()) {
				double next = iter.next().doubleValue();
				double diff = next - curr;
				rtrn.add(Double.valueOf(diff));
				curr = next;
			}
			
			// "Inter event time" between last and first event is time to first event + time after last event
			double timeToFirst = eventTimes.first().doubleValue();
			double timeAfterLast = siteExp.getTotalTime() - eventTimes.last().doubleValue();
			rtrn.add(Double.valueOf(timeToFirst + timeAfterLast));
			return rtrn;
			
		}
		
		/*
		 * Get the shortest time between two consecutive events
		 * Shortest time may be time to first event + time after last event
		 */
		@SuppressWarnings("unused")
		public double shortestInterEventTime() {
			return interEventTimes().first().doubleValue();
		}
		
		/*
		 * Get the longest time between two consecutive events
		 * Longest time may be time to first event + time after last event
		 */
		@SuppressWarnings("unused")
		public double longestInterEventTime() {
			return interEventTimes().last().doubleValue();
		}
		
		/*
		 * Get the collection of events at the site
		 */
		public TreeSet<Event> getEvents() {
			return siteEvents;
		}
		
		/*
		 * X coordinate of centroid
		 */
		public double getXCoord() {
			return xCoord;
		}
		
		/*
		 * Y coordinate of centroid
		 */
		public double getYCoord() {
			return yCoord;
		}
		
		/*
		 * Get the times of all events in the window
		 */
		public TreeSet<Double> getEventTimes(double windowMin, double windowMax, boolean inclusive) {
			TreeSet<Double> rtrn = new TreeSet<Double>();
			for(Event event : siteEvents) {
				double time = event.getRelTime();
				if(time < windowMin) continue;
				if(time > windowMax) continue;
				if(!inclusive) {
					if(time == windowMin) continue;
					if(time == windowMax) continue;
				}
				rtrn.add(Double.valueOf(time));
			}
			return rtrn;
		}
		
		/*
		 * Get all the event times as a tab delimited string
		 */
		@SuppressWarnings("unused")
		public String getEventTimesString() {
			TreeSet<Double> times = getEventTimes();
			String rtrn = "";
			for(Double time : times) {rtrn += time.toString() + "\t";}
			return rtrn;
		}
		
		/*
		 * Count events after the event within the time interval
		 */
		public int countEventsAfter(Event event) {
			double eventTime = event.getRelTime();
			int rtrn = 0;
			for(Event other : getEvents()) {
				if(other.equals(event)) continue;
				double otherTime = other.getRelTime();
				if(otherTime < eventTime) continue;
				if(otherTime - eventTime <= timeInterval) rtrn++;
			}
			return rtrn;
		}
		
		/*
		 * Get the next event or null if this is the last event
		 */
		public Event nextEvent(Event first) {
			Iterator<Event> iter = siteEvents.iterator();
			while(iter.hasNext()) {
				Event event = iter.next();
				if(event.equals(first)) {
					if(iter.hasNext()) return iter.next();
					return null;
				}
			}
			throw new IllegalArgumentException("Site does not have event.");
		}
		
		/*
		 * Whether there is another event after the event
		 */
		public boolean hasNextEvent(Event first) {
			return nextEvent(first) != null;
		}
		
		/*
		 * Get the time interval to the next event
		 */
		public double timeToNextEvent(Event first) {
			if(hasNextEvent(first)) {
				return nextEvent(first).getRelTime() - first.getRelTime();
			}
			throw new IllegalArgumentException("No next event");
		}
		
		/*
		 * Count the number of events in the window
		 */
		public int countEventsInWindow(double timeStart, double timeStop) {
			int rtrn = 0;
			for(Event event : getEvents()) {
				double time = event.getRelTime();
				if(time < timeStart) continue;
				if(time > timeStop) continue;
				rtrn++;
			}
			return rtrn;
		}
		
		/*
		 * Get all clusters of more than one event where each pair of consecutive events in the cluster happens within the time limit of each other
		 */
		public Collection<EventCluster> findEventClustersByInterEventTime() {
			Collection<EventCluster> rtrn = new ArrayList<EventCluster>();
			EventCluster currCluster = new EventCluster();
			for(Event event : siteEvents) {
				if(currCluster.isEmpty()) {
					currCluster.addEvent(event);
					continue;
				}
				double time = event.getRelTime();
				if(time - currCluster.lastTime() <= clusterMaxInterEventTime) {
					currCluster.addEvent(event);
					continue;
				}
				if(currCluster.numEvents() > 1) {
					rtrn.add(currCluster.copy());
				}
				currCluster.clear();
				currCluster.addEvent(event);
			}
			if(currCluster.numEvents() > 1) {
				rtrn.add(currCluster.copy());
			}
			return rtrn;
		}
		
		/*
		 * Get all possible clusters of consecutive events with no criteria
		 */
		@SuppressWarnings("unused")
		public Collection<EventCluster> getAllPossibleEventClusters() {
			Collection<EventCluster> rtrn = new ArrayList<EventCluster>();
			for(Event first : siteEvents) {
				Collection<Event> growingCluster = new TreeSet<Event>();
				for(Event second : siteEvents) {
					if(second.compareTo(first) < 0) continue;
					growingCluster.add(second);
					EventCluster tmpCluster = new EventCluster(growingCluster);
					rtrn.add(tmpCluster);
				}
			}
			return rtrn;
		}
				
		/*
		 * Get all possible clusters of consecutive events of the specified size
		 */
		public Collection<EventCluster> getAllEventClusters(int clusterSize) {
			Collection<EventCluster> rtrn = new ArrayList<EventCluster>();
			for(Event first : siteEvents) {
				Collection<Event> growingCluster = new TreeSet<Event>();
				for(Event second : siteEvents) {
					if(second.compareTo(first) < 0) continue;
					growingCluster.add(second);
					EventCluster tmpCluster = new EventCluster(growingCluster);
					if(tmpCluster.numEvents() == clusterSize) {
						rtrn.add(tmpCluster);
						break;
					}
				}
			}
			return rtrn;
		}

		
		/*
		 * Get the size of the largest cluster defined by inter event time
		 */
		public int largestClusterSizeByInterEventTime() {
			int rtrn = 0;
			Collection<EventCluster> clusters = findEventClustersByInterEventTime();
			for(EventCluster cluster : clusters) {
				if(cluster.numEvents() > rtrn) {
					rtrn = cluster.numEvents();
				}
			}
			return rtrn;
		}
		
		/*
		 * The maximal score of all possible clusters (consecutive events) of the specified size
		 * Calls the score method in EventCluster
		 */
		public double maxClusterScore(int clusterSize) {
			Collection<EventCluster> allClusters = getAllEventClusters(clusterSize);
			double rtrn = 0;
			for(EventCluster cluster : allClusters) {
				double score = cluster.score();
				if(score > rtrn) rtrn = score;
			}
			return rtrn;
		}
		
		/*
		 * Get all clusters such that the probability that any site in a random experiment
		 * has greater max cluster score is < alpha
		 * Filter clusters having max inter event time greater than cutoff
		 * Specify the cutoffs in parameter
		 */
		private ArrayList<EventCluster> getUncollapsedSignificantEventClusters(Map<Integer, Double> sigScoreCutoffs) {
			ArrayList<EventCluster> rtrn = new ArrayList<EventCluster>();
			for(int size = 2; size <= numEvents(); size++) {
				double scoreCutoff = sigScoreCutoffs.get(Integer.valueOf(size)).doubleValue();
				Collection<EventCluster> clusters = getAllEventClusters(size);
				for(EventCluster cluster : clusters) {
					if(cluster.score() > scoreCutoff && cluster.maxInterEventTime() <= clusterMaxInterEventTime) {
						rtrn.add(cluster);
					}
				}
			}
			return rtrn;
		}
		
		/*
		 * Get all clusters such that the probability that any site in a random experiment
		 * has greater max cluster score is < alpha
		 * Filter clusters having max inter event time greater than cutoff
		 */
		private ArrayList<EventCluster> getUncollapsedSignificantEventClusters() {
			if(!hasScoreCutoffs) {
				calculateScoreCutoffs();
			}
			return getUncollapsedSignificantEventClusters(scoreCutoffs);
		}

		
		/*
		 * Get all clusters such that the probability that a random site
		 * has greater max cluster score is < alpha
		 * Filter clusters having max inter event time greater than cutoff
		 * Then remove any cluster that is fully nested inside a larger significant cluster
		 */
		public ArrayList<EventCluster> significantEventClusters() {
			return collapseFullyNestedClusters(getUncollapsedSignificantEventClusters());
		}
		
		/*
		 * Get all clusters such that the probability that a random site
		 * has greater max cluster score is < alpha
		 * Filter clusters having max inter event time greater than cutoff
		 * Then remove any cluster that is fully nested inside a larger significant cluster
		 * Specify the cutoffs in parameter
		 */
		public ArrayList<EventCluster> significantEventClusters(Map<Integer, Double> sigScoreCutoffs) {
			return collapseFullyNestedClusters(getUncollapsedSignificantEventClusters(sigScoreCutoffs));
		}

		
		/*
		 * For each cluster size, get the 1 - alpha quantile of max cluster score in randomized sites
		 */
		private Map<Integer, Double> clusterScoreCutoffSingleSite() {
			Map<Integer, Double> rtrn = new TreeMap<Integer, Double>();
			for(int size = 2; size <= numEvents(); size++) {
				double[] randScores = new double[numRandomPermutations];
				for(int i=0; i < numRandomPermutations; i++) {
					Site rand = randomizeEventTimes();
					double randScore = rand.maxClusterScore(size);
					randScores[i] = randScore;
				}			
				Arrays.sort(randScores);
				rtrn.put(Integer.valueOf(size), Double.valueOf(Statistics.quantile(randScores, 1 - alpha)));
			}
			return rtrn;
		}
		
		/*
		 * Calculate significant score cutoff for cluster score by cluster size
		 */
		private void calculateScoreCutoffs() {
			System.out.println("Calculating site wide significant cluster score cutoffs for site " + siteID + "...");
			scoreCutoffs = clusterScoreCutoffSingleSite();
			hasScoreCutoffs = true;
		}
		
		/*
		 * Get significant score cutoff for cluster score by cluster size
		 */		
		@SuppressWarnings("unused")
		public Map<Integer, Double> scoreCutoffs() {
			if(hasScoreCutoffs) return scoreCutoffs;
			calculateScoreCutoffs();
			return scoreCutoffs;
		}
		
		/*
		 * Get a new site with same experiment and same number of events, with event times randomized
		 */
		public Site randomizeEventTimes() {
			double expTime = siteExp.getTotalTime();
			Site rtrn = new Site(siteExp, siteID + "_randomized_times", xCoord, yCoord);
			for(int i = 0; i < numEvents(); i++) {
				double randTime = Math.random() * expTime;
				rtrn.addEvent(randTime);
			}
			return rtrn;
		}
		
		/*
		 * (non-Javadoc)
		 * @see java.lang.Object#equals(java.lang.Object)
		 * Equal iff site ID, parent experiment, and number of events are equal
		 */
		@Override
		public boolean equals(Object o) {
			Site s = (Site)o;
			if(!s.getSiteId().equals(getSiteId())) return false;
			if(!s.getExperiment().equals(getExperiment())) return false;
			if(s.numEvents() != numEvents()) return false;
			return true;
		}
		
		/*
		 * (non-Javadoc)
		 * @see java.lang.Comparable#compareTo(java.lang.Object)
		 * Compare parent experiment, site ID, then number of events
		 */
		@Override
		public int compareTo(Site other) {
			if(this.equals(other)) return 0;
			if(!other.getExperiment().equals(getExperiment())) {
				return getExperiment().compareTo(other.getExperiment());
			}
			if(!other.getSiteId().equals(getSiteId())) {
				return getSiteId().compareTo(other.getSiteId());
			}
			if(other.numEvents() != numEvents()) {
				return Double.valueOf(numEvents()).compareTo(Double.valueOf(other.numEvents()));
			}
			throw new IllegalArgumentException("Can't compare sites with same experiment, same ID and same number of events:\n" + toString() + "\n" + other.toString());
		}

		/*
		 * (non-Javadoc)
		 * @see java.lang.Object#hashCode()
		 * Hash code of toString() representation
		 */
		@Override
		public int hashCode() {
			return toString().hashCode();
		}
		
		/*
		 * A cluster of events at the site
		 */
		private class EventCluster {
			
			/*
			 * The events in the cluster
			 */
			private TreeSet<Event> clusterEvents;
			
			/*
			 * Instantiate with no events
			 */
			public EventCluster() {
				clusterEvents = new TreeSet<Event>();
			}
			
			/*
			 * Construct with a collection of events
			 */
			public EventCluster(Collection<Event> events) {
				clusterEvents = new TreeSet<Event>();
				clusterEvents.addAll(events);
			}
			
			/*
			 * Clear the set of events
			 */
			public void clear() {
				clusterEvents.clear();
			}
			
			/*
			 * Get the events in the cluster
			 */
			public TreeSet<Event> getClusterEvents() {
				return clusterEvents;
			}
			
			/*
			 * Get the time of the last event in the cluster
			 */
			public double lastTime() {
				Iterator<Event> iter = clusterEvents.descendingIterator();
				return iter.next().getRelTime();
			}
			
			/*
			 * Get maximal time interval between two consecutive events
			 */
			public double maxInterEventTime() {
				if(numEvents() == 0) {
					throw new IllegalStateException("No events in cluster");
				}
				if(numEvents() == 1) return totalTimeSpan();
				double rtrn = 0;
				Iterator<Event> iter = clusterEvents.iterator();
				double currTime = iter.next().getRelTime();
				while(iter.hasNext()) {
					double nextTime = iter.next().getRelTime();
					double interEventTime = nextTime - currTime;
					if(interEventTime > rtrn) rtrn = interEventTime;
					currTime = nextTime;
				}
				return rtrn;
			}
			
			/*
			 * Get minimal time interval between two consecutive events
			 */
			@SuppressWarnings("unused")
			public double minInterEventTime() {
				if(numEvents() == 0) {
					throw new IllegalStateException("No events in cluster");
				}
				if(numEvents() == 1) return totalTimeSpan();
				double rtrn = totalTimeSpan();
				Iterator<Event> iter = clusterEvents.iterator();
				double currTime = iter.next().getRelTime();
				while(iter.hasNext()) {
					double nextTime = iter.next().getRelTime();
					double interEventTime = nextTime - currTime;
					if(interEventTime < rtrn) rtrn = interEventTime;
					currTime = nextTime;
				}
				return rtrn;
			}

			
			/*
			 * Get time from first event to last event
			 */
			public double totalTimeSpan() {
				return lastTime() - firstTime();
			}
			
			/*
			 * Get the time of the first event in the cluster
			 */
			public double firstTime() {
				Iterator<Event> iter = clusterEvents.iterator();
				return iter.next().getRelTime();
			}

			/*
			 * Reciprocal of the total time span of the cluster
			 */
			public double score() {
				if(numEvents() <= 1) return 0;
				return 1 / totalTimeSpan();
			}
			
			/*
			 * Get a copy of the cluster
			 */
			public EventCluster copy() {
				return new EventCluster(clusterEvents);
			}
			
			
			/*
			 * Whether every event in this cluster is in other cluster
			 */
			public boolean isNestedInside(EventCluster other) {
				TreeSet<Event> otherEvents = other.getClusterEvents();
				for(Event event : clusterEvents) {
					if(!otherEvents.contains(event)) return false;
				}
				return true;
			}
			
			/*
			 * Whether there are any events in the cluster
			 */
			public boolean isEmpty() {
				return clusterEvents.isEmpty();
			}
			
			/*
			 * Add an event to the cluster
			 */
			public void addEvent(Event event) {
				clusterEvents.add(event);
			}
			
			/*
			 * Get the number of events in the cluster
			 */
			public int numEvents() {
				return clusterEvents.size();
			}
			
			
			/*
			 * (non-Javadoc)
			 * @see java.lang.Object#toString()
			 * String consisting of string representations of all cluster events in order
			 */
			@Override
			public String toString() {
				String rtrn = siteID + "\t";
				for(Event event : clusterEvents) {
					rtrn += event.getRelTime() + "\t";
				}
				return rtrn;
			}
			
		}

	}
	
	/*
	 * Various types of scores for the site
	 */
	private class SiteScore {
		
		private double avgTimePerEvent;
		private double siteClusterPval;
		private double siteClusterScore;
		private Site site;
		
		public SiteScore(Site s) {
			site = s;
			calculateAvgTimePerEvent();
			calculateSiteClusterPval();
			calculateSiteClusterScore();
		}
		
		/*
		 * The average time per event at the site
		 */
		public double getAvgTimePerEvent() {
			return avgTimePerEvent;
		}
		
		/*
		 * The clusteriness score for the site
		 */
		public double getSiteClusterScore() {
			return siteClusterScore;
		}
		
		/*
		 * The p value for clusteriness of the site
		 */
		public double getSiteClusterPval() {
			return siteClusterPval;
		}
		
		private int NUM_BINS_OBSERVED_FREQUENCY_DIST = 10;
		
		/*
		 * P value of Pearson's chi square goodness of fit test
		 * Fitting observed frequency distribution of time between events to exponential distribution
		 */
		private void calculateSiteClusterPval() {

			ExponentialDistribution expDist = new ExponentialDistribution(avgTimePerEvent);
			TreeSet<Double> interEventTimes = site.interEventTimes();
			EmpiricalDistribution frequencyDist = new EmpiricalDistribution(interEventTimes, NUM_BINS_OBSERVED_FREQUENCY_DIST);
			
			// Number of frequencies in observed frequency distribution
			double numFrequencies = frequencyDist.numBins();
			
			// 1 + number of parameters in exponential distribution
			double reductionInDegreesOfFreedom = 2;
			
			// Degrees of freedom of test statistic
			double degreesOfFreedom = numFrequencies - reductionInDegreesOfFreedom;
			
			// Calculate test statistic
			double testStatistic = 0;
			for(int binNumber=0; binNumber < numFrequencies; binNumber++) {
				double observed = frequencyDist.getHistogram(binNumber);
				double binStart = frequencyDist.getBinStart(binNumber);
				double binEnd = frequencyDist.getBinEnd(binNumber);
				double expected = expDist.cumulativeProbability(binEnd) - expDist.cumulativeProbability(binStart);
				testStatistic += (Math.pow(observed - expected, 2)) / (expected);
			}
			
			// Get P value
			ChiSquaredDistribution chiSquaredDist = new ChiSquaredDistribution(degreesOfFreedom);
			siteClusterPval = 1 - chiSquaredDist.cumulativeProbability(testStatistic);
			
		}
		
		/*
		 * Score is log10(site cluster p value)
		 */
		private void calculateSiteClusterScore() {
			siteClusterScore = Math.log10(siteClusterPval);
		}
		
		
		private void calculateAvgTimePerEvent() {
			avgTimePerEvent = site.getExperiment().getTotalTime() / site.numEvents();
		}
		
	}
	

	
	/*
	 * A release event at a site
	 */
	private class Event implements Comparable<Event> {
		
		/*
		 * The parent experiment
		 */
		private Experiment eventExp;
		
		/*
		 * The parent site
		 */
		private Site eventSite;
		
		/*
		 * The time of the event
		 */
		private double eventTime;
		
		/*
		 * Construct with parent experiment, parent site and event time
		 */
		public Event(Experiment experiment, Site site, double relTime) {
			eventExp = experiment;
			eventSite = site;
			eventTime = relTime;
		}

		/*
		 * Get the parent experiment
		 */
		public Experiment getExperiment() {
			return eventExp;
		}

		/*
		 * Get the parent site
		 */
		public Site getSite() {
			return eventSite;
		}

		/*
		 * Get the event time
		 */
		public double getRelTime() {
			return eventTime;
		}

		/*
		 * (non-Javadoc)
		 * @see java.lang.Object#equals(java.lang.Object)
		 * Equal iff parent experiment, parent site and event time are equal
		 */
		@Override
		public boolean equals(Object o) {
			Event e = (Event)o;
			if(!e.getExperiment().equals(getExperiment())) return false;
			if(!e.getSite().equals(getSite())) return false;
			return e.getRelTime() == getRelTime();
		}

		/*
		 * (non-Javadoc)
		 * @see java.lang.Comparable#compareTo(java.lang.Object)
		 * Compare parent experiment, then parent site, then event time
		 */
		@Override
		public int compareTo(Event other) {
			if(!getExperiment().equals(other.getExperiment())) {
				return getExperiment().compareTo(other.getExperiment());
			}
			if(!getSite().equals(other.getSite())) {
				return getSite().compareTo(other.getSite());
			}
			return Double.valueOf(getRelTime()).compareTo(Double.valueOf(other.getRelTime()));
		}
		
		/*
		 * (non-Javadoc)
		 * @see java.lang.Object#toString()
		 * String with helpful information
		 */
		@Override
		public String toString() {
			return "experiment_name:" + eventExp.getName() + ";site_name:" + eventSite.getSiteId() + ";rel_time:" + Double.valueOf(eventTime).toString();
		}

		/*
		 * (non-Javadoc)
		 * @see java.lang.Object#hashCode()
		 * Hash code of toString() representation
		 */
		@Override
		public int hashCode() {
			return toString().hashCode();
		}

		
	}
	
	/*
	 * All the sites within a certain radius of a point
	 */
	private class Neighborhood {
		
		private Collection<Site> sites;
		
		/*
		 * A neighborhood comprising a certain radius around a site
		 */
		@SuppressWarnings("unused")
		public Neighborhood(Experiment exp, Site centerSite, Distance dist, double radius) {
			
			if(!exp.getSites().contains(centerSite)) {
				throw new IllegalArgumentException("Center site must belong to the experiment.");
			}
			sites = getAllSitesWithinRadius(exp, centerSite.getXCoord(), centerSite.getYCoord(), radius, dist);
			
		}
		
		/*
		 * Construct neighborhood containing the sites by taking centroid and radius of collection and including all sites from experiment that lie within radius of centroid
		 */
		@SuppressWarnings("unused")
		public Neighborhood(Experiment exp, Collection<Site> siteCollection, Distance dist) {
			double radius = Dataset.this.radius(siteCollection);
			double centroidX = Dataset.centroidXcoord(siteCollection);
			double centroidY = Dataset.centroidYcoord(siteCollection);
			sites = Dataset.getAllSitesWithinRadius(exp, centroidX, centroidY, radius, dist);
		}
		
		/*
		 * A neighborhood comprising a certain radius around a point
		 */
		@SuppressWarnings("unused")
		public Neighborhood(Experiment exp, double centerX, double centerY, Distance dist, double radius) {
			
			sites = getAllSitesWithinRadius(exp, centerX, centerY, radius, dist);
			
		}

		/*
		 * The total number of events at sites in the neighborhood
		 */
		public int numEvents() {
			int rtrn = 0;
			for(Site site : sites) {
				rtrn += site.numEvents();
			}
			return rtrn;
		}
		
		/*
		 * The total number of events in the time interval at sites in the neighborhood
		 */
		@SuppressWarnings("unused")
		public int numEvents(double timeStart, double timeStop) {
			int rtrn = 0;
			for(Site site : sites) {
				rtrn += site.getEventTimes(timeStart, timeStop, true).size();
			}
			return rtrn;
		}
		
		/*
		 * The number of sites
		 */
		public int numSites() {
			return sites.size();
		}
		
		/*
		 * The average number of events per site
		 */
		public double meanEventsPerSite() {
			return (double) numEvents() / (double) numSites();
		}
		
		/*
		 * The total number of events at all sites divided by the total time measured at all sites
		 */
		public double meanTimePerEvent() {
			int numEvents = numEvents();
			double expTime = getExperiment().getTotalTime();
			double totalTime = numSites() * expTime;
			return numEvents / totalTime;
		}
		
		/*
		 * Get all the sites in neighborhoon
		 */
		public Collection<Site> getSites() {
			return sites;
		}
		
		/*
		 * Get the experiment
		 */
		public Experiment getExperiment() {
			return sites.iterator().next().getExperiment();
		}
		
		/*
		 * Distance from centroid to furthest site
		 */
		public double radius() {
			return Dataset.this.radius(sites);
		}
		
		/*
		 * X coordinate of centroid of all sites
		 */
		public double centroidXcoord() {
			return Dataset.centroidXcoord(sites);
		}
		
		/*
		 * Y coordinate of centroid of all sites
		 */
		public double centroidYcoord() {
			return Dataset.centroidXcoord(sites);
		}
		
		
		/*
		 * Whether all the sites in other are in this neighborhood
		 */
		@SuppressWarnings("unused")
		private boolean fullyContains(Neighborhood other) {
			for(Site otherSite : other.getSites()) {
				if(!sites.contains(otherSite)) return false;
			}
			return true;
		}
		
		/*
		 * (non-Javadoc)
		 * @see java.lang.Object#equals(java.lang.Object)
		 * 
		 */
		@Override
		public boolean equals(Object o) {
			Neighborhood n = (Neighborhood)o;
			return sites.equals(n.getSites());
		}

		/*
		 * A string consisting of tab separated names of the fields contained in the toString() representation
		 */
		@SuppressWarnings("unused")
		public String getToStringFieldNamesAsString() {
			String rtrn = "Experiment\t";
			rtrn += "Num_sites\t";
			rtrn += "Centroid_X\t";
			rtrn += "Centroid_Y\t";
			rtrn += "Radius\t";
			rtrn += "Mean_events_per_site\t";
			rtrn += "Mean_time_per_event\t";
			rtrn += "Site_IDs";
			return rtrn;
		}
		
		/*
		 * (non-Javadoc)
		 * @see java.lang.Object#toString()
		 * String with helpful information
		 */
		@Override
		public String toString() {
			String rtrn = getExperiment().getName() + "\t";
			rtrn += numSites() + "\t";
			rtrn += centroidXcoord() + "\t";
			rtrn += centroidYcoord() + "\t";
			rtrn += radius() + "\t";
			rtrn += meanEventsPerSite() + "\t";
			rtrn += meanTimePerEvent() + "\t";
			for(Site site : sites) {
				rtrn += site.getSiteId() + ";";
			}
			return rtrn;
		}

		/*
		 * (non-Javadoc)
		 * @see java.lang.Object#hashCode()
		 * Hash code of toString() representation
		 */
		@Override
		public int hashCode() {
			return toString().hashCode();
		}
		
	}
	
	/*
	 * A distance between two points
	 */
	private interface Distance {
		
		/*
		 * Distance between two sites
		 */
		public double getDistance(Site site1, Site site2);
		
		/*
		 * Distance between two points described by (x,y) coordinates
		 */
		public double getDistance(double xCoord1, double yCoord1, double xCoord2, double yCoord2);
		
		/*
		 * Distance between a point described by (x,y) coordinates and a site
		 */
		public double getDistance(double xCoord, double yCoord, Site site);
	}
	
	/*
	 * A Euclidian distance
	 */
	private class EuclidianDistance implements Distance {

		public EuclidianDistance() {}
		
		@Override
		public double getDistance(Site site1, Site site2) {
			double x1 = site1.getXCoord();
			double y1 = site1.getYCoord();
			double x2 = site2.getXCoord();
			double y2 = site2.getYCoord();
			return Math.sqrt(Math.pow(x1-x2,2) + Math.pow(y1-y2,2));
		}
		
		@Override
		public double getDistance(double xCoord1, double yCoord1, double xCoord2, double yCoord2) {
			return Math.sqrt(Math.pow(xCoord1-xCoord2,2) + Math.pow(yCoord1-yCoord2,2));
		}

		@Override
		public double getDistance(double xCoord, double yCoord, Site site) {
			double xCoord1 = site.getXCoord();
			double yCoord1 = site.getYCoord();
			return getDistance(xCoord, yCoord, xCoord1, yCoord1);
		}
		
	}
	

}
