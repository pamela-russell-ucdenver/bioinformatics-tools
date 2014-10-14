package fasta;

import general.CommandLineParser;

import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

import org.apache.log4j.Logger;

import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;

/**
 * @author prussell
 *
 */
public class FastaSizes {
	
	private static Logger logger = Logger.getLogger(FastaSizes.class.getName());
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-i", "Fasta file", true);
		p.addStringArg("-o", "Output table of sizes", true);
		p.parse(args);
		String input = p.getStringArg("-i");
		String output = p.getStringArg("-o");
		
		FastaSequenceIO fsio = new FastaSequenceIO(input);
		List<Sequence> seqs = fsio.loadAll();
		
		logger.info("Writing sizes to file " + output + "...");
		
		FileWriter w = new FileWriter(output);
		for(Sequence seq : seqs) {
			w.write(seq.getId() + "\t" + seq.getLength() + "\n");
		}
		w.close();

		logger.info("All done.");
		
	}

}
