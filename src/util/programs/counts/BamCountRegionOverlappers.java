package util.programs.counts;

import java.io.IOException;

import org.apache.log4j.Logger;

import net.sf.samtools.util.CloseableIterator;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.MappedFragment;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.AbstractAnnotationCollection;
import guttmanlab.core.annotationcollection.BAMFragmentCollectionFactory;
import guttmanlab.core.annotationcollection.FeatureCollection;

/**
 * Count the number of bam records that overlap features in an annotation
 * @author prussell
 *
 */
public class BamCountRegionOverlappers {
	
	private FeatureCollection<Gene> features;
	private AbstractAnnotationCollection<? extends MappedFragment> data;
	private static Logger logger = Logger.getLogger(BamCountRegionOverlappers.class.getName());
	
	/**
	 * @param bamFile Bam file
	 * @param geneBed Bed file of features to check for overlappers
	 * @param chrSizes Chromosome size file
	 * @throws IOException
	 */
	public BamCountRegionOverlappers(String bamFile, String geneBed, String chrSizes) throws IOException {
		features = (FeatureCollection<Gene>) BEDFileIO.loadFromFile(geneBed, chrSizes);
		data = BAMFragmentCollectionFactory.createFromBam(bamFile);
	}
	
	/**
	 * @return Total number of bam records that overlap a feature in the annotation
	 */
	public int getTotalOverlappers() {
		CloseableIterator<? extends MappedFragment> iter = data.sortedIterator();
		int numDone = 0;
		int rtrn = 0;
		while(iter.hasNext()) {
			numDone++;
			if(features.overlaps(iter.next())) {
				rtrn++;
			}
			if(numDone % 10000000 == 0) {
				logger.info("Finished " + numDone + " records of which " + rtrn + " overlap features.");
			}
		}
		iter.close();
		return rtrn;
	}
	
}
