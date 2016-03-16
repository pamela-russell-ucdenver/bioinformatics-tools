package util.programs.bed;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;

import org.apache.log4j.Logger;

import util.FileUtils;
import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.FeatureCollection;
import guttmanlab.core.util.CommandLineParser;

public class MultipleBedFileOverlapper {
	
	
	private Map<String, Map<Integer, Integer>> overlappersByPosition;
	private boolean countIntrons;
	private static Logger logger = Logger.getLogger(MultipleBedFileOverlapper.class.getName());
	
	/**
	 * @param bedFileList File containing list of bed files to overlap
	 * @param chrSizeFile Chromosome size file
	 * @throws IOException
	 */
	private MultipleBedFileOverlapper(String bedFileList, String chrSizeFile) throws IOException {
		this(bedFileList, false, chrSizeFile);
	}
	
	/**
	 * @param bedFileList File containing list of bed files to overlap
	 * @param includeIntrons Include introns so annotation is considered to be entire gene span
	 * @param chrSizeFile Chromosome size file
	 * @throws IOException
	 */
	private MultipleBedFileOverlapper(String bedFileList, boolean includeIntrons, String chrSizeFile) throws IOException {
		countIntrons = includeIntrons;
		initializeOverlapperCountsFromBedFiles(FileUtils.fileLinesAsList(bedFileList), chrSizeFile);
	}
	
	private void initializeOverlapperCountsFromBedFiles(Collection<String> bedFiles, String chrSizeFile) throws IOException {
		logger.info("");
		overlappersByPosition = new TreeMap<String, Map<Integer, Integer>>();
		int totalPos = 0;
		for(String bedFile : bedFiles) {
			logger.info("Reading file " + bedFile + "...");
			Map<String, FeatureCollection<Gene>> genes = BEDFileIO.loadFromFileByReferenceName(bedFile, chrSizeFile);
			for(String chr : genes.keySet()) {
				if(!overlappersByPosition.containsKey(chr)) {
					overlappersByPosition.put(chr, new TreeMap<Integer, Integer>());
					logger.info("Added chromosome " + chr + ".");
				}
				for(Gene gene : genes.get(chr)) {
					Gene geneToUse = countIntrons ? new Gene(new SingleInterval(gene.getReferenceName(), gene.getReferenceStartPosition(), gene.getReferenceEndPosition(), gene.getOrientation())) : gene;
					Iterator<SingleInterval> iter = geneToUse.getBlocks();
					while(iter.hasNext()) {
						SingleInterval block = iter.next();
						for(int i = block.getReferenceStartPosition(); i < block.getReferenceEndPosition(); i++) {
							Integer pos = Integer.valueOf(i);
							if(!overlappersByPosition.get(chr).containsKey(pos)) {
								overlappersByPosition.get(chr).put(pos, Integer.valueOf(0));
							}
							overlappersByPosition.get(chr).put(pos, Integer.valueOf(overlappersByPosition.get(chr).get(pos).intValue() + 1));
							totalPos++;
						}
					}
				}
			}
			logger.info("Done reading file " + bedFile + ". Total coverage of all positions: " + totalPos + ".");
		}
	}
	
	private Collection<SingleInterval> getContiguousIntervalsWithMinCount(int minCount) {
		logger.info("");
		logger.info("Getting contiguous intervals with min count " + minCount + "...");
		Collection<SingleInterval> rtrn = new ArrayList<SingleInterval>();
		for(String chr : overlappersByPosition.keySet()) {
			int numIntervals = 0;
			boolean withinInterval = false;
			int currentIntervalStart = 0;
			int currentIntervalMax = 0;
			Iterator<Integer> iter = overlappersByPosition.get(chr).keySet().iterator();
			while(iter.hasNext()) {
				int pos = iter.next().intValue();
				int count = overlappersByPosition.get(chr).get(Integer.valueOf(pos)).intValue();
				if(count >= minCount) {
					if(withinInterval) {
						if(pos == currentIntervalMax + 1) {
							// Continue with current interval
							currentIntervalMax = pos;
							continue;
						} else {
							// Write previous interval and start new interval here
							rtrn.add(new SingleInterval(chr, currentIntervalStart, currentIntervalMax));
							numIntervals++;
							currentIntervalStart = pos;
							currentIntervalMax = pos;
						}
					} else {
						// Start new interval here
						currentIntervalStart = pos;
						currentIntervalMax = pos;
						withinInterval = true;
					}
				} else {
					if(withinInterval) {
						// Write previous interval and exit interval
						rtrn.add(new SingleInterval(chr, currentIntervalStart, currentIntervalMax));
						numIntervals++;
						withinInterval = false;
					} else {
						continue;
					}
				}
			}
			logger.info("Got " + numIntervals + " intervals on chromosome " + chr + ".");
		}
		return rtrn;
	}
	
	/**
	 * @param intervals Intervals to write
	 * @param minIntervalSize Min size to write an interval
	 * @param bothOrientations Write two copies of the interval, one in each orientation
	 * @param outBed Output bed file
	 * @throws IOException
	 */
	private static void writeToFile(Collection<SingleInterval> intervals, int minIntervalSize, boolean bothOrientations, String outBed) throws IOException {
		logger.info("");
		logger.info("Writing to file " + outBed + "...");
		FileWriter w = new FileWriter(outBed);
		for(SingleInterval interval : intervals) {
			if(interval.size() < minIntervalSize) {
				continue;
			}
			String chr = interval.getReferenceName();
			int start = interval.getReferenceStartPosition();
			int end = interval.getReferenceEndPosition();
			if(bothOrientations) {
				String plusName = chr + ":" + start + "-" + end + ":+";
				String minusName = chr + ":" + start + "-" + end + ":-";
				SingleInterval plus = new SingleInterval(chr, start, end, Strand.POSITIVE, plusName);
				SingleInterval minus = new SingleInterval(chr, start, end, Strand.NEGATIVE, minusName);
				w.write(plus.toBED() + "\n");
				w.write(minus.toBED() + "\n");
			} else {
				String name = chr + ":" + start + "-" + end;
				SingleInterval si = new SingleInterval(chr, start, end, Strand.POSITIVE, name);
				w.write(si.toBED() + "\n");
			}
		}
		w.close();
	}
	
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-l", "File containing list of bed files to overlap", true);
		p.addBooleanArg("-s", "Count full gene span including introns", false, false);
		p.addBooleanArg("-b", "Write intervals in both orientations", false, false);
		p.addIntArg("-mo", "Min number of overlappers per position", true);
		p.addIntArg("-ms", "Min size of interval to write", false, 0);
		p.addStringArg("-o", "Output bed file", true);
		p.addStringArg("-c", "Chromosome size file", true);
		p.parse(args);
		String listFile = p.getStringArg("-l");
		boolean fullSpan = p.getBooleanArg("-s");
		boolean bothOrientations = p.getBooleanArg("-b");
		int minOverlappers = p.getIntArg("-mo");
		int minSize = p.getIntArg("-ms");
		String outBed = p.getStringArg("-o");
		String chrSizeFile = p.getStringArg("-c");
		
		MultipleBedFileOverlapper m = new MultipleBedFileOverlapper(listFile, fullSpan, chrSizeFile);
		Collection<SingleInterval> allIntervals = m.getContiguousIntervalsWithMinCount(minOverlappers);
		writeToFile(allIntervals, minSize, bothOrientations, outBed);
		
		logger.info("");
		logger.info("All done.");
	}

}
