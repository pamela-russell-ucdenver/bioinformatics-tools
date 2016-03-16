package util;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.FeatureCollection;
import guttmanlab.core.coordinatespace.CoordinateSpace;

import java.io.IOException;
import java.util.Iterator;
import java.util.Map;

import org.apache.log4j.Logger;

import net.sf.samtools.util.CloseableIterator;

public class OverlapUtils {
	
	private Map<String, FeatureCollection<Gene>> genesByChr;
	private CoordinateSpace coordSpace;
	@SuppressWarnings("unused")
	private static Logger logger = Logger.getLogger(OverlapUtils.class.getName());
	
	/**
	 * @param bedFile Gene annotation in bed format
	 * @param referenceSizes Reference genome size file. Line format: chr   size
	 * @throws IOException
	 */
	public OverlapUtils(String bedFile, String referenceSizes) throws IOException {
		loadGenes(bedFile, referenceSizes);
		coordSpace = new CoordinateSpace(referenceSizes);
	}
	
	private void loadGenes(String bedFile, String referenceSizes) throws IOException {
		genesByChr = BEDFileIO.loadFromFileByReferenceName(bedFile, referenceSizes);
	}
	
	/**
	 * Get genes overlapping an interval
	 * @param chr Interval chr
	 * @param start Interval start
	 * @param end Interval end
	 * @return Genes overlapping the interval
	 */
	public FeatureCollection<Gene> getOverlappers(String chr, int start, int end) {
		SingleInterval interval = new SingleInterval(chr, start, end);
		FeatureCollection<Gene> rtrn = new FeatureCollection<Gene>(coordSpace);
		if(genesByChr.get(chr).overlaps(interval)) {
			CloseableIterator<Gene> overlappers = genesByChr.get(chr).sortedIterator(interval, false);
			while(overlappers.hasNext()) {
				Gene gene = overlappers.next();
				if(gene.overlaps(interval)) { // extra check for gene blocks
					rtrn.add(gene);
				}
			}
			overlappers.close();
		}
		return rtrn;
	}

	/**
	 * Get genes overlapping a position
	 * @param chr Position chr
	 * @param pos Position coord
	 * @return Genes overlapping the position
	 */
	public FeatureCollection<Gene> getOverlappers(String chr, int pos) {
		return getOverlappers(chr, pos, pos+1);
	}

	
	/**
	 * Get exon number overlapping a position
	 * @param annot The blocked feature
	 * @param chr Position chr
	 * @param genomicPos Position
	 * @return Overlapping exon number (counting left to right in genomic coordinates)
	 */
	public static int exonNumber(Annotation annot, String chr, int genomicPos) {
		if(!chr.equals(annot.getReferenceName())) {
			throw new IllegalArgumentException("Position not on same chromosome as annotation");
		}
		Iterator<SingleInterval> iter = annot.getBlocks();
		int exonNum = -1;
		while(iter.hasNext()) {
			SingleInterval interval = iter.next();
			exonNum++;
			if(genomicPos >= interval.getReferenceStartPosition() && genomicPos < interval.getReferenceEndPosition()) {
				return exonNum;
			}
		}
		throw new IllegalArgumentException("Position doesn't overlap an exon");
	}
	
	
	/**
	 * Get exon number overlapping a position, counting from 5' end of the feature
	 * @param annot The blocked feature
	 * @param chr Position chr
	 * @param genomicPos Position
	 * @return Overlapping exon number (counting from 5' end of feature)
	 */
	public static int exonNumberFrom5Prime(Annotation annot, String chr, int genomicPos) {
		int e = exonNumber(annot, chr, genomicPos);
		if(annot.getOrientation().equals(Strand.POSITIVE)) {
			return e;
		}
		if(annot.getOrientation().equals(Strand.NEGATIVE)) {
			return annot.getNumberOfBlocks() - e - 1;
		}
		throw new IllegalArgumentException("Strand must be + or -");
	}
	
	
	

}