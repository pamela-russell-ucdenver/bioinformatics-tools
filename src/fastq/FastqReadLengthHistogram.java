package fastq;

import guttmanlab.core.util.CommandLineParser;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import org.apache.log4j.Logger;

import broad.core.math.EmpiricalDistribution;

public class FastqReadLengthHistogram {
	
	private static Logger logger = Logger.getLogger(FastqReadLengthHistogram.class.getName());
	
	private static EmpiricalDistribution getDistribution(String fastqFile, int maxBin) throws IOException {
		EmpiricalDistribution dist = new EmpiricalDistribution(maxBin, 0, maxBin, true);
		BufferedReader reader = new BufferedReader(new FileReader(fastqFile));
		int numRead = 0;
		while(reader.ready()) {
			String line = reader.readLine();
			numRead++;
			if(numRead % 4 == 2) {
				dist.add(line.length());
			}
			if(numRead % 400000 == 0) {
				int n = numRead / 4;
				logger.info("Finished " + n + " reads.");
			}
		}
		reader.close();
		return dist;
	}
	
	private static void writeHistogram(EmpiricalDistribution dist, String outFile) throws IOException {
		FileWriter w = new FileWriter(outFile);
		int numBins = dist.numBins();
		for(int i = 0; i < numBins; i++) {
			double binStart = dist.getBinStart(i);
			double binEnd = dist.getBinEnd(i);
			double binCount = dist.getHistogram(i);
			w.write(binStart + "-" + binEnd + "\t" + binCount + "\n");
		}
		w.close();
	}
	
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-f", "fastq file", true);
		p.addIntArg("-m", "max bin", true);
		p.addStringArg("-o", "out histogram file", true);
		p.parse(args);
		String fastq = p.getStringArg("-f");
		int maxBin = p.getIntArg("-m");
		String outFile = p.getStringArg("-o");
		
		writeHistogram(getDistribution(fastq, maxBin), outFile);
		
		logger.info("");
		logger.info("All done.");
		
	}

}
