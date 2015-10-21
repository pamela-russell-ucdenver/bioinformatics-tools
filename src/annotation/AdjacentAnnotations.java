package annotation;

import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.AnnotationCollection;
import guttmanlab.core.annotationcollection.FeatureCollection;
import guttmanlab.core.coordinatespace.CoordinateSpace;
import guttmanlab.core.util.CommandLineParser;

import java.io.IOException;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.log4j.Logger;

import net.sf.samtools.util.CloseableIterator;

/**
 * Look in a target annotation set for features adjacent to features in a query set
 * @author prussell
 *
 */
public class AdjacentAnnotations {
	
	private static Logger logger = Logger.getLogger(AdjacentAnnotations.class.getName());
	
	/**
	 * Load targets as tree sets by chromosome
	 * @param bed Bed file of targets
	 * @param refSizes Reference size file
	 * @return Map of chromosome to sorted set of targets on chromosome
	 * @throws IOException
	 */
	private static Map<String, TreeSet<Gene>> loadAsTreeSets(String bed, String refSizes) throws IOException {
		logger.info("");
		logger.info("Loading targets as sorted sets...");
		Map<String, FeatureCollection<Gene>> features = BEDFileIO.loadFromFileByReferenceName(bed, refSizes);
		Map<String, TreeSet<Gene>> rtrn = new TreeMap<String, TreeSet<Gene>>();
		for(String chr : features.keySet()) {
			logger.info(chr);
			rtrn.put(chr, new TreeSet<Gene>());
			CloseableIterator<Gene> iter = features.get(chr).sortedIterator();
			while(iter.hasNext()) {
				rtrn.get(chr).add(iter.next());
			}
			iter.close();
		}
		return rtrn;
	}
	
	/**
	 * Get greatest element in target set strictly less than query gene
	 * @param query Query gene
	 * @param targets Target set
	 * @return Greatest element in target set strictly less than query gene
	 */
	private static Gene previous(Gene query, Map<String, TreeSet<Gene>> targets) {
		return targets.get(query.getReferenceName()).lower(query);
	}

	/**
	 * Get least element in target set strictly greater than query gene
	 * @param query Query gene
	 * @param targets Target set
	 * @return Least element in target set strictly greater than query gene
	 */
	private static Gene next(Gene query, Map<String, TreeSet<Gene>> targets) {
		return targets.get(query.getReferenceName()).higher(query);
	}
	
	/**
	 * Get all previous and next elements for a collection of queries
	 * @param queries Query genes
	 * @param targets Sorted target genes by chromosome
	 * @param refSizes Reference size file
	 * @return FeatureCollection of all previous and next elements for a collection of queries
	 */
	private static FeatureCollection<Gene> getAllAdjacent(AnnotationCollection<Gene> queries, Map<String, TreeSet<Gene>> targets, String refSizes) {
		FeatureCollection<Gene> rtrn = new FeatureCollection<Gene>(new CoordinateSpace(refSizes));
		CloseableIterator<Gene> iter = queries.sortedIterator();
		while(iter.hasNext()) {
			Gene gene = iter.next();
			Gene prev = previous(gene, targets);
			if(prev == null) {
				logger.info("No previous gene for " + gene.getName());
			} else {
				rtrn.add(prev);
			}
			Gene next = next(gene, targets);
			if(next == null) {
				logger.info("No next gene for " + gene.getName());
			} else {
				rtrn.add(next);
			}
		}
		iter.close();
		return rtrn;
	}
	
	/**
	 * Get all previous and next elements for a collection of queries and write to file
	 * @param queryBed Bed file of query genes
	 * @param targetBed Bed file of target genes
	 * @param outputBed Output bed file of adjacent genes
	 * @param refSizes Reference size file
	 * @throws IOException
	 */
	private static void writeAllAdjacent(String queryBed, String targetBed, String outputBed, String refSizes) throws IOException {
		logger.info("");
		logger.info("Getting adjacent genes...");
		AnnotationCollection<Gene> queries = BEDFileIO.loadFromFile(queryBed, refSizes);
		Map<String, TreeSet<Gene>> targets = loadAsTreeSets(targetBed, refSizes);
		FeatureCollection<Gene> adjacent = getAllAdjacent(queries, targets, refSizes);
		logger.info("");
		logger.info("Writing to bed file...");
		BEDFileIO.writeToFile(adjacent, outputBed);
	}
	
	
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-q", "Query bed file", true);
		p.addStringArg("-t", "Target bed file", true);
		p.addStringArg("-o", "Output bed file", true);
		p.addStringArg("-r", "Reference size file", true);
		p.parse(args);
		String queryBed = p.getStringArg("-q");
		String targetBed = p.getStringArg("-t");
		String outputBed = p.getStringArg("-o");
		String refSizes = p.getStringArg("-r");
		
		writeAllAdjacent(queryBed, targetBed, outputBed, refSizes);
		
		logger.info("");
		logger.info("All done.");
	}
	
	
	
}
