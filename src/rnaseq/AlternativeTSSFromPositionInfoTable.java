package rnaseq;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeSet;

import org.apache.log4j.Logger;

import broad.core.parser.CommandLineParser;
import broad.pda.annotation.BEDFileParser;
import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Gene;
import nextgen.core.annotation.Annotation.Strand;
import nextgen.core.annotation.BasicAnnotation;


public class AlternativeTSSFromPositionInfoTable extends MaxFragmentEndPileupFromPositionInfoTable {
	
	private static Logger logger = Logger.getLogger(AlternativeTSSFromPositionInfoTable.class.getName());
	
	private AlternativeTSSFromPositionInfoTable() {
		super();
	}
	
	/**
	 * Get all positions passing a count threshold
	 * Positions within a distance cutoff of already included position are not included
	 * Recursively eliminate buffer zone around included position
	 * If no positions passing count cutoff, return empty collection
	 * Also enforce MIN_PILEUP_SINGLE_POS
	 * @param region Region
	 * @param tableFile File of counts by position
	 * @param colNum Column number to get counts from
	 * @param countCutoff Min count to include pileup
	 * @param distCutoff Min distance between returned positions
	 * @return All positions with pileup passing count cutoff, separated by the buffer
	 * @throws IOException
	 */
	private Collection<PositionAndCount> getPositionsPassingCutoffs(Annotation region, String tableFile, int colNum, float countCutoff, int distCutoff) throws IOException {
		Collection<PositionAndCount> rtrn = new TreeSet<PositionAndCount>();
		PositionAndCount max = getMaxPosition(region, tableFile, colNum, false);
		if(max.getCount() < countCutoff) return rtrn;
		if(max.getCount() < MIN_PILEUP_SINGLE_POS) return rtrn;
		rtrn.add(max);
		int maxPos = max.getPos();
		Annotation buffer = new BasicAnnotation(region.getChr(), maxPos - distCutoff, maxPos + distCutoff, region.getOrientation());
		Annotation newRegion = region.minus(buffer);
		if(newRegion.getSize() < 1) return rtrn;
		if(newRegion.equals(region)) return rtrn;
		rtrn.addAll(getPositionsPassingCutoffs(newRegion, tableFile, colNum, countCutoff, distCutoff));
		return rtrn;
	}
	
	/**
	 * Get all positions within a % cutoff of the max pileup in the region
	 * If the max is less than MIN_PILEUP_SINGLE_POS, returns null
	 * Also, positions within a distance cutoff of already included position are not included
	 * Recursively eliminate buffer zone around included position
	 * Only return pileups passing MIN_PILEUP_SINGLE_POS
	 * @param region Region
	 * @param tableFile File of counts by position
	 * @param colNum Column number to get counts from
	 * @param pctCutoffRelToMax Min % relative to max pileup
	 * @param distCutoff Min distance between returned positions
	 * @return Max pileup position and all others within the cutoff of the max, separated by the buffer
	 * @throws IOException
	 */
	private Collection<PositionAndCount> getTopPostionsRelToMax(Annotation region, String tableFile, int colNum, float pctCutoffRelToMax, int distCutoff) throws IOException {
		if(!(pctCutoffRelToMax > 0 && pctCutoffRelToMax < 1)) {
			throw new IllegalArgumentException("Cutoff compared to max must be between 0 and 1");
		}
		PositionAndCount max = getMaxPosition(region, tableFile, colNum, false);
		int maxCount = max.getCount();
		if(maxCount < MIN_PILEUP_SINGLE_POS) {
			logger.warn("No pileups over " + MIN_PILEUP_SINGLE_POS + " for region " + region.getName());
			return null;
		}
		float countCutoff = pctCutoffRelToMax * maxCount;
		return getPositionsPassingCutoffs(region, tableFile, colNum, countCutoff, distCutoff);
	}

	private static Annotation tssFlank(Annotation gene, int upstream, int downstream) {
		String rtrnName = gene.getName() + "_tss_flank_" + upstream + "_" + downstream;
		String rtrnChr = gene.getReferenceName();
		Strand rtrnStrand = gene.getOrientation();
		int rtrnStart = -1;
		int rtrnEnd = -1;
		switch(rtrnStrand) {
		case NEGATIVE:
			rtrnStart = gene.getEnd() - downstream;
			rtrnEnd = gene.getEnd() + upstream;
			return new BasicAnnotation(rtrnChr, rtrnStart, rtrnEnd, rtrnStrand, rtrnName);
		case POSITIVE:
			rtrnStart = gene.getStart() - upstream;
			rtrnEnd = gene.getStart() + downstream;
			return new BasicAnnotation(rtrnChr, rtrnStart, rtrnEnd, rtrnStrand, rtrnName);
		case UNKNOWN:
			throw new IllegalArgumentException("Strand must be positive or negative");
		default:
			throw new IllegalArgumentException("Strand must be positive or negative");
		}
	}
	
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-b", "Bed annotation", true);
		p.addStringArg("-t", "Position info table", true);
		p.addStringArg("-o", "Output bed file of alternative TSS's", true);
		p.addIntArg("-u", "Distance upstream of TSS to look for alternative start sites", true);
		p.addIntArg("-d", "Distance downstream of TSS to look for alternative start sites", true);
		p.addFloatArg("-pc", "Pct cutoff relative to max pileup", true);
		p.addIntArg("-dc", "Min distance between returned alternative starts", true);
		p.addIntArg("-c", "Zero-based column number in table for 5' reads starting at position", true);
		p.parse(args);
		String bedFile = p.getStringArg("-b");
		String table = p.getStringArg("-t");
		String output = p.getStringArg("-o");
		int upstream = p.getIntArg("-u");
		int downstream = p.getIntArg("-d");
		float pctCutoffRelToMax = p.getFloatArg("-pc");
		int distCutoff = p.getIntArg("-dc");
		int startersCol5p = p.getIntArg("-c");
		
		Map<String, Collection<Gene>> genes = BEDFileParser.loadDataByChr(bedFile);
		
		AlternativeTSSFromPositionInfoTable a = new AlternativeTSSFromPositionInfoTable();
		
		FileWriter w = new FileWriter(output);
		for(String chr : genes.keySet()) {
			logger.info(chr);
			for(Gene gene : genes.get(chr)) {
				try {
					Collection<PositionAndCount> pileups = a.getTopPostionsRelToMax(tssFlank(gene, upstream, downstream), 
							table, startersCol5p, pctCutoffRelToMax, distCutoff);
					if(pileups == null) continue;
					for(PositionAndCount pileup : pileups) {
						String line = gene.getChr() + "\t";
						line += pileup.getPos() + "\t";
						line += Integer.valueOf(pileup.getPos() + 1).toString() + "\t";
						line += gene.getName() + "_alt_tss_" + pileup.getPos() + "\t";
						line += pileup.getCount() + "\t";
						line += gene.getOrientation().toString();
						w.write(line + "\n");
					}
				} catch (IllegalStateException e) {
					e.printStackTrace();
					continue;
				}
			}
		}
		w.close();
		
		logger.info("");
		logger.info("All done.");
		
	}

	
	
}
