/**
 * 
 */
package test;

import guttmanlab.core.util.CommandLineParser;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;

import org.apache.log4j.Logger;

import broad.pda.annotation.BEDFileParser;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.BasicAnnotation;
import nextgen.core.annotation.Gene;
import nextgen.core.coordinatesystem.GenomicSpace;
import nextgen.core.coordinatesystem.TranscriptomeSpace;
import nextgen.core.model.AlignmentModel;
import nextgen.core.model.TranscriptomeSpaceAlignmentModel;
import nextgen.core.model.score.ScanStatisticScore;
import nextgen.core.readFilters.GenomicSpanFilter;



/**
 * @author prussell
 *
 */
public class TestFragmentCounts {

	private TranscriptomeSpaceAlignmentModel transcriptomeData;
	private AlignmentModel genomicData;
	private Map<String, Collection<Gene>> genes;
	private static int DEFAULT_MAX_GENOMIC_SPAN = 300000;
	private static Logger logger = Logger.getLogger(TestFragmentCounts.class.getName());

	private TestFragmentCounts(String bamFile, String bedFile, String chrSizeFile) throws IOException {
		genes = BEDFileParser.loadDataByChr(new File(bedFile));
		transcriptomeData = new TranscriptomeSpaceAlignmentModel(bamFile, new TranscriptomeSpace(genes));
		transcriptomeData.addFilter(new GenomicSpanFilter(DEFAULT_MAX_GENOMIC_SPAN));
		//transcriptomeData.addFilter(new NumHitsFilter(1));
		genomicData = new AlignmentModel(bamFile, new GenomicSpace(chrSizeFile), false);
		genomicData.addFilter(new GenomicSpanFilter(DEFAULT_MAX_GENOMIC_SPAN));
		//genomicData.addFilter(new NumHitsFilter(1));
	}
	
	@SuppressWarnings("unused")
	private ScanStatisticScore scoreWindow(Annotation window) {
		ScanStatisticScore score = new ScanStatisticScore(transcriptomeData, window, false);
		double regionLength = window.getSize();
		double regionTotal = transcriptomeData.getCount(window);
		score.setRegionLength(regionLength);
		score.setRegionTotal(regionTotal);
		score.refreshScanPvalue(transcriptomeData);
		return score;
	}

	private void printCounts(String chr, int start, int end) {
		Annotation window = new BasicAnnotation(chr, start, end);
		logger.info(window.toBED());
		//ScanStatisticScore score = scoreWindow(window);
		//double scoreCount = score.getCount();
		//logger.info("scan_statistic_score_count\t" + scoreCount);
		//double scoreRegionTotal = score.getRegionTotal();
		//logger.info("scan_statistic_score_region_total\t" + scoreRegionTotal);
		double transcriptomeCount = transcriptomeData.getCount(window);
		logger.info("transcriptome_count\t" + transcriptomeCount);
		double genomeCount = genomicData.getCount(window);
		logger.info("genome_count\t" + genomeCount);		
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-b", "Bam file", true);
		p.addStringArg("-a", "Bed annotation", true);
		p.addStringArg("-cs", "Chromosome size file", true);
		p.addStringArg("-c", "Chromosome name", true);
		p.addIntArg("-s", "Start", true);
		p.addIntArg("-e", "End", true);
		p.parse(args);
		
		String bamFile = p.getStringArg("-b");
		String bedFile = p.getStringArg("-a");
		String chrSizeFile = p.getStringArg("-cs");
		TestFragmentCounts t = new TestFragmentCounts(bamFile,bedFile,chrSizeFile);
		
		String chr = p.getStringArg("-c");
		int start = p.getIntArg("-s");
		int end = p.getIntArg("-e");
		
		t.printCounts(chr, start, end);

	}

}
