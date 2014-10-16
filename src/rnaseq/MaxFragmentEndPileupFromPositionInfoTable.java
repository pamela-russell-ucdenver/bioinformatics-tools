package rnaseq;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

import org.apache.log4j.Logger;

import broad.core.parser.CommandLineParser;
import broad.core.parser.StringParser;
import broad.pda.annotation.BEDFileParser;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Annotation.Strand;
import nextgen.core.annotation.Gene;

public class MaxFragmentEndPileupFromPositionInfoTable {
	
	private StringParser stringParser;
	private static Logger logger = Logger.getLogger(MaxFragmentEndPileupFromPositionInfoTable.class.getName());
	
	private MaxFragmentEndPileupFromPositionInfoTable() {
		stringParser = new StringParser();
	}
	
	private class PositionAndCount {
		
		private String chr;
		private int pos;
		private int count;
		
		public PositionAndCount(String chrName, int position, int num) {
			chr = chrName;
			pos = position;
			count = num;
		}
		
		@SuppressWarnings("unused")
		public String getChr() {
			return chr;
		}
		
		public int getPos() {
			return pos;
		}
		
		@SuppressWarnings("unused")
		public int getCount() {
			return count;
		}
		
	}
	
	
	private int getCount(String tableLine, int colNum) {
		stringParser.parse(tableLine);
		return stringParser.asInt(colNum);
	}
	
	private int getPos(String tableLine) {
		stringParser.parse(tableLine);
		return stringParser.asInt(1);
	}

	private String getChr(String tableLine) {
		stringParser.parse(tableLine);
		return stringParser.asString(0);
	}
	
	private boolean matchesChr(Annotation region, String tableLine) {
		String chr = region.getChr();
		return getChr(tableLine).equals(chr);
	}
	
	private boolean withinSpan(Annotation region, String tableLine) {
		int start = region.getStart();
		int end = region.getEnd();
		int pos = getPos(tableLine);
		return pos >= start && pos < end;
	}
	
	/**
	 * Move the reader to the second line of the annotation and return the first line as a string
	 * @param tableReader Buffered reader for position info table
	 * @param region Region of interest
	 * @return First line of region
	 * @throws IOException 
	 */
	private String seekToStartPos(BufferedReader tableReader, Annotation region) throws IOException {
		while(tableReader.ready()) {
			String line = tableReader.readLine();
			if(!matchesChr(region, line)) {
				continue;
			}
			if(withinSpan(region, line)) {
				return line;
			}
		}
		throw new IllegalStateException("Never found span of " + region.getName());
	}
	
	private Map<Integer, Integer> getCounts(Annotation region, String tableFile, int colNum) throws IOException {
		Map<Integer, Integer> rtrn = new TreeMap<Integer, Integer>();
		BufferedReader reader = new BufferedReader(new FileReader(tableFile));
		String firstLine = seekToStartPos(reader, region);
		rtrn.put(Integer.valueOf(getPos(firstLine)), Integer.valueOf(getCount(firstLine, colNum)));
		while(reader.ready()) {
			String line = reader.readLine();
			if(!withinSpan(region, line)) {
				return rtrn;
			}
			rtrn.put(Integer.valueOf(getPos(line)), Integer.valueOf(getCount(line, colNum)));
		}
		return rtrn;
	}
	
	private PositionAndCount getMaxPosition(Annotation region, String tableFile, int colNum) throws IOException {
		Map<Integer, Integer> counts = getCounts(region, tableFile, colNum);
		int maxPos = Integer.MIN_VALUE;
		int maxCount = Integer.MIN_VALUE;
		for(Integer pos : counts.keySet()) {
			int val = counts.get(pos).intValue();
			if(val > maxCount) {
				maxPos = pos.intValue();
				maxCount = val;
			}
		}
		return new PositionAndCount(region.getChr(), maxPos, maxCount);
	}
	
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-b", "Bed annotation", true);
		p.addStringArg("-t", "Position info table", true);
		p.addStringArg("-o", "Output bed", true);
		p.parse(args);
		String bedFile = p.getStringArg("-b");
		String table = p.getStringArg("-t");
		String output = p.getStringArg("-o");
		
		Map<String, Collection<Gene>> genes = BEDFileParser.loadDataByChr(bedFile);
		
		int startersCol5p = 4;
		int startersCol3p = 5;
		MaxFragmentEndPileupFromPositionInfoTable m = new MaxFragmentEndPileupFromPositionInfoTable();
		
		FileWriter w = new FileWriter(output);
		for(String chr : genes.keySet()) {
			logger.info(chr);
			for(Gene gene : genes.get(chr)) {
				try {
					PositionAndCount p5p = m.getMaxPosition(gene, table, startersCol5p);
					PositionAndCount p3p = m.getMaxPosition(gene, table, startersCol3p);
					
					int end5p = p5p.getPos();
					int end3p = p3p.getPos();
					Strand strand = gene.getStrand();
					if(strand.equals(Strand.UNKNOWN)) {
						throw new IllegalArgumentException("Strand must be known");
					}
					int startPos = strand.equals(Strand.POSITIVE) ? end5p : end3p;
					int endPos = strand.equals(Strand.POSITIVE) ? end3p : end5p;
					if(startPos >= endPos) {
						throw new IllegalStateException("Start is greater than end (" + startPos + " > " + endPos + ")");
					}
					
					String line = gene.getChr() + "\t";
					line += startPos + "\t";
					line += endPos + "\t";
					line += gene.getName() + "\t";
					line += "0.0\t";
					line += strand.toString();

					w.write(line + "\n");
				} catch(IllegalStateException e) {
					logger.warn("SKIPPING GENE " + gene.getName() + " " + e.getMessage());
				}
			}
		}
		w.close();
		
		logger.info("");
		logger.info("All done.");
		
	}

}
