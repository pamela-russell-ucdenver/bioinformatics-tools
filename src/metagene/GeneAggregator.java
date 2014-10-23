package metagene;

import guttmanlab.core.util.CommandLineParser;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.apache.log4j.Logger;

import broad.pda.annotation.BEDFileParser;

import nextgen.core.annotation.Gene;

/**
 * @author prussell
 *
 */
public class GeneAggregator {
	
	private GenePartition genePartition;
	private RegionDataType regionDataType;
	private Map<String, Collection<Gene>> genes;
	private Map<GeneComponent, Collection<Gene>> regionsByComponent;
	private Map<GeneComponent, List<Double>> summedDataByComponent;
	private Map<GeneComponent, Map<Gene, List<Double>>> data;
	
	private static int DEFAULT_NORMALIZED_REGION_SIZE = 100;
	private static double DEFAULT_NORMALIZED_VECTOR_SUM = 1;
	private static Logger logger = Logger.getLogger(GeneAggregator.class.getName());
	
	/**
	 * @param components Components of gene partition
	 * @param dataType Data type
	 * @param bedAnnotation Bed annotation of gene set
	 * @throws IOException
	 */
	public GeneAggregator(Collection<GeneComponent> components, RegionDataType dataType, String bedAnnotation) throws IOException {
		this(new GenePartition(components), dataType, BEDFileParser.loadDataByChr(new File(bedAnnotation)));
	}
	
	/**
	 * @param partition Gene partition
	 * @param dataType Data type
	 * @param genesByChr Gene set
	 */
	public GeneAggregator(GenePartition partition, RegionDataType dataType, Map<String, Collection<Gene>> genesByChr) {
		genePartition = partition;
		regionDataType = dataType;
		genes = genesByChr;
		regionsByComponent = genePartition.getPartition(genes);
		summedDataByComponent = new HashMap<GeneComponent, List<Double>>();
		data = new HashMap<GeneComponent, Map<Gene, List<Double>>>();
		for(GeneComponent component : genePartition.getComponents()) {
			data.put(component, new TreeMap<Gene, List<Double>>());
		}
		
	}
	
	/**
	 * For each gene component, write table of data vector for each region
	 * @param outfilePrefix Prefix for output files
	 * @throws IOException
	 */
	@SuppressWarnings("unused")
	private void writeTables(String outfilePrefix) throws IOException {
		writeTables(outfilePrefix, false, 0);
	}
	
	/**
	 * For each gene component, write table of (possibly normalized) data vector for each region
	 * @param outfilePrefix Prefix for output files
	 * @param normalizeOverallSum Whether to normalize to a particular sum for each region
	 * @param newSum The new sum to normalize to for each component
	 * @throws IOException
	 */
	private void writeTables(String outfilePrefix, boolean normalizeOverallSum, double newSum) throws IOException {
		logger.info("");
		logger.info("Writing tables for each region...");
		for(GeneComponent component : genePartition.getComponents()) {
			String outTable = outfilePrefix + "_" + component.getName();
			Collection<Gene> regions = regionsByComponent.get(component);
			logger.info("Writing table of " + regions.size() + " regions for gene component " + component.getName() + " to file " + outTable + "...");
			FileWriter w = new FileWriter(outTable);
			int numDone = 0;
			for(Gene region : regions) {
				numDone++;
				if(numDone % 1000 == 0) {
					logger.info("Finished " + numDone + " regions.");
				}
				String line = region.getName() + "\t";
				List<Double> dataList = getData(region, component);
				if(normalizeOverallSum) {
					List<Double> normalizedData = Util.normalizeSum(dataList, newSum);
					line += Util.listToString(normalizedData);
				} else {
					line += Util.listToString(dataList);
				}
				w.write(line + "\n");
			}
			w.close();
		}
		logger.info("Done writing tables.");
	}
	
	/**
	 * Get the vector of data values for the region and component
	 * Get from cache or cache after computing
	 * @param region The region
	 * @param component Gene component
	 * @return Data values for region
	 * @throws IOException
	 */
	private List<Double> getData(Gene region, GeneComponent component) throws IOException {
		if(data.get(component).containsKey(region)) {
			return data.get(component).get(region);
		}
		List<Double> rtrn = regionDataType.getData(region, component.reverseDataIfMinusOrientation());
		data.get(component).put(region, rtrn);
		return rtrn;
	}
	
	/**
	 * Write overall total data vector for one gene component to file
	 * @param component Gene component
	 * @param outFile Output file
	 * @throws IOException
	 */
	@SuppressWarnings("unused")
	private void writeTotalOverAllRegions(GeneComponent component, String outFile) throws IOException {
		writeTotalOverAllRegions(component, outFile, false, 0);
	}
	
	/**
	 * Write (possibly normalized) overall total data vector for one gene component to file
	 * @param component Gene component
	 * @param outFile Output file
	 * @param normalizeOverallSum Whether to normalize to a particular overall sum for each component
	 * @param newSum The new sum to normalize to for each component
	 * @throws IOException
	 */
	private void writeTotalOverAllRegions(GeneComponent component, String outFile, boolean normalizeOverallSum, double newSum) throws IOException {
		logger.info("Writing total over all regions for gene component " + component.getName() + " to file " + outFile + "...");
		FileWriter w = new FileWriter(outFile);
		if(normalizeOverallSum) {
			logger.info("Normalizing vector sum to " + newSum);
			List<Double> normalizedTotal = getNormalizedTotalOverAllRegions(component, newSum);
			w.write(Util.listToTableByPos(normalizedTotal));
		} else {
			List<Double> dataTotal = getTotalOverAllRegions(component);
			w.write(Util.listToTableByPos(dataTotal));
		}
		w.close();
		logger.info("Done writing file.");
	}
	
	/**
	 * Write overall total data vector for each gene component to file
	 * @param outFilePrefix Output file prefix
	 * @throws IOException
	 */
	@SuppressWarnings("unused")
	private void writeTotalOverAllRegions(String outFilePrefix) throws IOException {
		writeTotalOverAllRegions(outFilePrefix, false, 0);
	}
	
	/**
	 * Write (possibly normalized) overall total data vector for each gene component to file
	 * @param outFilePrefix Output file prefix
	 * @param normalizeOverallSum Whether to normalize to a particular overall sum for each component
	 * @param newSum The new sum to normalize to for each component
	 * @throws IOException
	 */
	private void writeTotalOverAllRegions(String outFilePrefix, boolean normalizeOverallSum, double newSum) throws IOException {
		logger.info("");
		logger.info("Writing total over all regions for each gene component...");
		for(GeneComponent component : genePartition.getComponents()) {
			String outFile = outFilePrefix + "_" + component.getName();
			writeTotalOverAllRegions(component, outFile, normalizeOverallSum, newSum);
		}
		logger.info("Done writing totals.");
	}
	
	/**
	 * Get normalized total of data over all regions in the component
	 * @param component The gene component to use
	 * @param newOverallTotal The overall sum of the normalized result
	 * @return A normalized list representing the sum of data for all regions
	 * @throws IOException
	 */
	private List<Double> getNormalizedTotalOverAllRegions(GeneComponent component, double newOverallTotal) throws IOException {
		logger.info("Getting normalized total over all regions for component " + component.getName());
		return Util.normalizeSum(getTotalOverAllRegions(component), newOverallTotal);
	}
	
	/**
	 * Get total of data over all regions in the component
	 * @param component The gene component to use
	 * @return A list representing the sum of data for all regions
	 * @throws IOException
	 */
	private List<Double> getTotalOverAllRegions(GeneComponent component) throws IOException {
		logger.info("Getting total over all regions for component " + component.getName());
		// Check if the result is already cached
		if(summedDataByComponent.containsKey(component)) {
			return summedDataByComponent.get(component);
		}
		Collection<Gene> regions = regionsByComponent.get(component);
		logger.info("");
		logger.info("Calculating data total over " + regions.size() + " regions in component " + component.getName() + "...");
		Collection<List<Double>> dataList = new ArrayList<List<Double>>();
		int numDone = 0;
		for(Gene region : regions) {
			numDone++;
			if(numDone % 1000 == 0) {
				logger.info("Finished " + numDone + " regions.");
			}
			dataList.add(getData(region, component));
		}
		List<Double> rtrn = Util.sum(dataList);
		// Cache the result
		summedDataByComponent.put(component, rtrn);
		return rtrn;
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-b", "Bam file", true);
		p.addStringArg("-g", "Gene annotation", true);
		p.addIntArg("-s", "Normalized size for each region", false, DEFAULT_NORMALIZED_REGION_SIZE);
		p.addStringArg("-ot", "Output file prefix for tables of all regions", false, null);
		p.addStringArg("-o", "Output file prefix for aggregated data", false, null);
		p.addBooleanArg("-n", "Normalize data vectors to a particular sum", false, false);
		p.addDoubleArg("-ns", "Normalized sum for each data vector (only works if -n is true", false, DEFAULT_NORMALIZED_VECTOR_SUM);
		p.addBooleanArg("--cds", "Include CDS", false, false);
		p.addBooleanArg("--utr3", "Include 3'-UTR", false, false);
		p.addBooleanArg("--utr5", "Include 5'-UTR", false, false);
		p.addBooleanArg("--gene-3-prime-end", "Include 3' end of transcript", false, false);
		p.addIntArg("-n3g", "Number of bases at 3' end of gene", false, 100);
		p.addBooleanArg("--downstream-of-splice-junction", "Include region downstream of splice junction", false, false);
		p.addIntArg("-nds", "Number of bases downstream of splice junction", false, 100);
		p.addBooleanArg("--upstream-of-splice-junction", "Include region upstream of splice junction", false, false);
		p.addIntArg("-nus", "Number of bases upstream of splice junction", false, 100);
		p.addBooleanArg("--5-prime-end-of-intron", "Include 5' region of intron", false, false);
		p.addIntArg("-n5i", "Number of bases at 5' end of intron", false, 100);
		p.addBooleanArg("--3-prime-end-of-intron", "Include 3' region of intron", false, false);
		p.addIntArg("-n3i", "Number of bases at 3' end of intron", false, 100);
		p.addBooleanArg("--introns", "Include introns", false, false);
		p.addBooleanArg("--whole-gene", "Include whole gene", false, false);
		p.addIntArg("-m5", "Minimum size of 5'-UTR to include", false, 0);
		p.addIntArg("-m3", "Minimum size of 3'-UTR to include", false, 0);
		
		
		p.parse(args);
		
		String bamFile = p.getStringArg("-b");
		String bedFile = p.getStringArg("-g");
		int regionSize = p.getIntArg("-s");
		String outTablePrefix = p.getStringArg("-ot");
		String outFilePrefix = p.getStringArg("-o");
		boolean normalizeSums = p.getBooleanArg("-n");
		double normalizedSum = p.getDoubleArg("-ns");
		boolean useCds = p.getBooleanArg("--cds");
		boolean useUtr3 = p.getBooleanArg("--utr3");
		boolean useUtr5 = p.getBooleanArg("--utr5");
		boolean useGene3PrimeEnd = p.getBooleanArg("--gene-3-prime-end");
		int numBasesGene3PrimeEnd = p.getIntArg("-n3g");
		boolean useDownstreamSpliceJunction = p.getBooleanArg("--downstream-of-splice-junction");
		int numBasesDownstreamSpliceJunction = p.getIntArg("-nds");
		boolean useUpstreamSpliceJunction = p.getBooleanArg("--upstream-of-splice-junction");
		int numBasesUpstreamSpliceJunction = p.getIntArg("-nus");
		boolean useIntron5primeEnd = p.getBooleanArg("--5-prime-end-of-intron");
		int numBasesIntron5primeEnd = p.getIntArg("-n5i");
		boolean useIntron3primeEnd = p.getBooleanArg("--3-prime-end-of-intron");
		int numBasesIntron3primeEnd = p.getIntArg("-n3i");
		boolean useIntrons = p.getBooleanArg("--introns");
		boolean useWholeGene = p.getBooleanArg("--whole-gene");
		int min5utr = p.getIntArg("-m5");
		int min3utr = p.getIntArg("-m3");
		
		
		CoverageByNormalizedPosition c = new CoverageByNormalizedPosition(bamFile, bedFile, regionSize);
		Collection<GeneComponent> components = new ArrayList<GeneComponent>();
		Collection<GeneComponent> requireBothUTRs = new ArrayList<GeneComponent>();
		Collection<GeneComponent> require3UTR = new ArrayList<GeneComponent>();
		Collection<GeneComponent> require5UTR = new ArrayList<GeneComponent>();
		GeneComponent3UTR require3utr = new GeneComponent3UTR(min3utr);
		GeneComponent5UTR require5utr = new GeneComponent5UTR(min5utr);
		requireBothUTRs.add(require3utr);
		requireBothUTRs.add(require5utr);
		require3UTR.add(require3utr);
		require5UTR.add(require5utr);
		if(useCds) components.add(new GeneComponentCDS(requireBothUTRs));
		if(useUtr3) components.add(new GeneComponent3UTR(require5UTR, min3utr));
		if(useUtr5) components.add(new GeneComponent5UTR(require3UTR, min5utr));
		if(useGene3PrimeEnd) components.add(new GeneComponent3PrimeEnd(numBasesGene3PrimeEnd));
		//if(useGene5PrimeEnd) components.add(new GeneComponent5PrimeEnd(numBasesGene5PrimeEnd));
		if(useDownstreamSpliceJunction) components.add(new GeneComponentDownstreamOfSpliceJunction(numBasesDownstreamSpliceJunction));
		if(useUpstreamSpliceJunction) components.add(new GeneComponentUpstreamOfSpliceJunction(numBasesUpstreamSpliceJunction));
		if(useIntron5primeEnd) components.add(new GeneComponentIntron5PrimeEnd(numBasesIntron5primeEnd));
		if(useIntron3primeEnd) components.add(new GeneComponentIntron3PrimeEnd(numBasesIntron3primeEnd));
		if(useIntrons) components.add(new GeneComponentIntrons());
		if(useWholeGene) components.add(new GeneComponentWholeGene());
		
		if(components.isEmpty()) {
			throw new IllegalArgumentException("Must choose at least one gene component to include.");
		}
		logger.info("");
		logger.info("Using gene components:");
		for(GeneComponent co : components) {
			logger.info(co.getName());
		}
		GeneAggregator a = new GeneAggregator(components, c, bedFile);
		
		if(outTablePrefix != null) {
			a.writeTables(outTablePrefix, normalizeSums, normalizedSum);
		}

		if(outFilePrefix != null) {
			a.writeTotalOverAllRegions(outFilePrefix, normalizeSums, normalizedSum);
		}

		
	}

}
