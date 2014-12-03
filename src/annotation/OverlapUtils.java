package annotation;

import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.FeatureCollection;
import guttmanlab.core.coordinatespace.CoordinateSpace;

import java.io.IOException;
import java.util.Map;

import net.sf.samtools.util.CloseableIterator;

public class OverlapUtils {
	
	private Map<String, FeatureCollection<Gene>> genesByChr;
	private CoordinateSpace coordSpace;

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
	 * Count number of positions in list that overlap some exon from annotation
	 * Write list of overlaps to a file
	 * @param positionList List file. Line format: chr   pos
	 * @param outFile Output file
	 * @throws IOException 
	 */
	public FeatureCollection<Gene> getOverlappers(String chr, int pos) {
		SingleInterval interval = new SingleInterval(chr, pos, pos + 1);
		FeatureCollection<Gene> rtrn = new FeatureCollection<Gene>(coordSpace);
		if(genesByChr.get(chr).overlaps(interval)) {
			CloseableIterator<Gene> overlappers = genesByChr.get(chr).sortedIterator(interval, false);
			while(overlappers.hasNext()) {
				Gene gene = overlappers.next();
				rtrn.add(gene);
			}
			overlappers.close();
		}
		return rtrn;
	}


}