package util.programs.text;

import guttmanlab.core.util.CommandLineParser;
import guttmanlab.core.util.StringParser;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeSet;

import org.apache.log4j.Logger;


/**
 * Make a histogram of number of occurrences of a string in a list of strings
 * @author prussell
 *
 */
public class StringFrequencyHistogram {
	
	private static Logger logger = Logger.getLogger(StringFrequencyHistogram.class.getName());
	
	private static void addString(Map<String, Integer> counts, String s) {
		if(counts.containsKey(s)) {
			counts.put(s, Integer.valueOf(counts.get(s).intValue() + 1));
		} else {
			counts.put(s, Integer.valueOf(1));
		}
	}
	
	@SuppressWarnings("unused")
	private void makeCounts(String file) throws IOException {
		makeCounts(file, -1);
	}
	
	private static Map<String, Integer> makeCounts(String file, int colNum) throws IOException {
		logger.info("Counting string frequencies...");
		Map<String, Integer> counts = new HashMap<String, Integer>();
		int numDone = 0;
		FileReader r = new FileReader(file);
		BufferedReader b = new BufferedReader(r);
		StringParser s = new StringParser();
		while(b.ready()) {
			if(numDone % 100000 == 0) {
				logger.info("Finished " + numDone + " lines.");
			}
			String line = b.readLine();
			if(colNum < 0) {
				addString(counts, line);
			} else {
				s.parse(line);
				String col = s.asString(colNum);
				addString(counts, col);
			}
			numDone++;
		}
		r.close();
		b.close();
		return counts;
	}
	
	private TreeSet<StringWithCount> sort(Map<String, Integer> counts) {
		logger.info("");
		logger.info("Sorting...");
		TreeSet<StringWithCount> rtrn = new TreeSet<StringWithCount>();
		for(String s : counts.keySet()) {
			rtrn.add(new StringWithCount(s, counts.get(s).intValue()));
		}
		return rtrn;
	}
	
	private static void write(Collection<StringWithCount> counts, String outFile) throws IOException {
		logger.info("");
		logger.info("Writing " + counts.size() + " counts to " + outFile);
		FileWriter w = new FileWriter(outFile);
		for(StringWithCount s : counts) {
			w.write(s.str + "\t" + s.cnt + "\n");
		}
		w.close();
	}
	
	private class StringWithCount implements Comparable<StringWithCount> {
		
		private String str;
		private int cnt;
		
		public StringWithCount(String string, int count) {
			str = string;
			cnt = count;
		}
		
		@Override
		public int compareTo(StringWithCount other) {
			if(cnt != other.cnt) {
				return other.cnt - cnt;
			}
			return other.str.compareTo(str);
		}
		
		
		
	}
	
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-i", "Input file or table", true);
		p.addStringArg("-o", "Output histogram table", true);
		p.addIntArg("-c", "Zero-based column number to count, or omit for whole line", false, -1);
		p.parse(args);
		String input = p.getStringArg("-i");
		String output = p.getStringArg("-o");
		int col = p.getIntArg("-c");
		
		write(new StringFrequencyHistogram().sort(makeCounts(input,col)), output);
		
		logger.info("");
		logger.info("All done.");
		
	}

}
