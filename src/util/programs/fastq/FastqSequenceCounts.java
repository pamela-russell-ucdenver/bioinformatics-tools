package util.programs.fastq;

import guttmanlab.core.util.CommandLineParser;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.apache.log4j.Logger;

public class FastqSequenceCounts {
	
	private Map<String, Integer> counts;
	private static Logger logger = Logger.getLogger(FastqSequenceCounts.class.getName());
	
	private FastqSequenceCounts() {
		counts = new HashMap<String, Integer>();
	}
	
	private void addSequence(String seq) {
		if(!counts.containsKey(seq)) {
			counts.put(seq, Integer.valueOf(1));
		} else {
			counts.put(seq, Integer.valueOf(counts.get(seq).intValue()+1));
		}
	}
	
	private void processFqFile(String fastqFile, int minReadLen) throws IOException {
		BufferedReader reader = new BufferedReader(new FileReader(fastqFile));
		int numRead = 0;
		while(reader.ready()) {
			String line = reader.readLine();
			numRead++;
			if(numRead % 4 == 2) {
				if(line.length() >= minReadLen) addSequence(line);
			}
			if(numRead % 400000 == 0) {
				int n = numRead / 4;
				logger.info("Finished " + n + " reads. There are " + counts.size() + " different sequences.");
			}
		}
		reader.close();
	}
	
	private void writeTable(String outFile) throws IOException {
		FileWriter w = new FileWriter(outFile);
		for(String seq : counts.keySet()) {
			w.write(seq + "\t" + counts.get(seq) + "\n");
		}
		w.close();
	}
	
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-i", "Fastq file", true);
		p.addStringArg("-o", "Output table of sequence counts", true);
		p.addIntArg("-ml", "Min read length to count", false, 0);
		p.parse(args);
		
		String fastq = p.getStringArg("-i");
		String out = p.getStringArg("-o");
		int minLen = p.getIntArg("-ml");
		
		FastqSequenceCounts counts = new FastqSequenceCounts();
		counts.processFqFile(fastq, minLen);
		counts.writeTable(out);
		
	}
	
}
