package metagene;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.log4j.Logger;

import broad.pda.annotation.BEDFileParser;

import nextgen.core.annotation.Annotation.Strand;
import nextgen.core.annotation.Gene;
import nextgen.core.coordinatesystem.TranscriptomeSpace;
import nextgen.core.model.TranscriptomeSpaceAlignmentModel;
import nextgen.core.readFilters.GenomicSpanFilter;

/**
 * @author prussell
 *
 */
public class CoverageByNormalizedPosition implements RegionDataType {

	private TranscriptomeSpaceAlignmentModel data;
	private int size;
	private static Logger logger = Logger.getLogger(CoverageByNormalizedPosition.class.getName());
	private static int DEFAULT_GENOMIC_SPAN_FILTER = 300000;
	
	/**
	 * @param bamFile Bam file
	 * @param transcriptomeBedFile Bed annotation for transcriptome space
	 * @param normalizedRegionSize Normalized size for each region
	 * @throws IOException
	 */
	public CoverageByNormalizedPosition(String bamFile, String transcriptomeBedFile, int normalizedRegionSize) throws IOException {
		this(new TranscriptomeSpaceAlignmentModel(bamFile, new TranscriptomeSpace(BEDFileParser.loadDataByChr(new File(transcriptomeBedFile)))), normalizedRegionSize);
	}
	
	/**
	 * @param alignmentData Alignment data
	 * @param normalizedRegionSize Size to normalize all regions to
	 */
	public CoverageByNormalizedPosition(TranscriptomeSpaceAlignmentModel alignmentData, int normalizedRegionSize) {
		logger.info("");
		logger.info("Instantiating CoverageByNormalizedPosition object...");
		data = alignmentData;
		data.addFilter(new GenomicSpanFilter(DEFAULT_GENOMIC_SPAN_FILTER));
		size = normalizedRegionSize;
	}
	
	@Override
	public List<Double> getData(Gene region, boolean reverseIfMinusOrientation) throws IOException {
		List<Double> fullCounts = data.getPositionCountList(region);
		if(reverseIfMinusOrientation && region.getOrientation().equals(Strand.NEGATIVE)) {
			ArrayList<Double> copy = new ArrayList<Double>();
			int s = fullCounts.size();
			for(int i=0; i<s; i++) {
				copy.add(new Double(fullCounts.get(i).doubleValue()));
			}
			for(int i=0; i<s; i++) {
				fullCounts.set(i, copy.get(s-i-1));
			}
		}
		return Util.expandOrContractList(fullCounts, size);
	}
	
	@Override
	public double getSummary(Gene region) {
		throw new UnsupportedOperationException("Not implemented");
	}

}
