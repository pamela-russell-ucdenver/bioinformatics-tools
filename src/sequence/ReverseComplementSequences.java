package sequence;

import general.CommandLineParser;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import org.apache.log4j.Logger;

import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;

/**
 * @author prussell
 *
 */
public class ReverseComplementSequences {
	
	private Collection<Sequence> sequences;
	private static Logger logger = Logger.getLogger(ReverseComplementSequences.class.getName());
	
	private ReverseComplementSequences(String inputFasta) throws IOException {
		FastaSequenceIO fsio = new FastaSequenceIO(inputFasta);
		sequences = fsio.loadAll();
	}
	
	private void writeReverseComplements(String outFasta, int lineLength) throws IOException {
		logger.info("Writing to file " + outFasta);
		List<Sequence> rcs = new ArrayList<Sequence>();
		for(Sequence seq : sequences) {
			Sequence rc = new Sequence(seq.getId() + "_RC");
			rc.setSequenceBases(seq.getSequenceBases());
			rc.reverse();
			rcs.add(rc);
		}
		FastaSequenceIO fsio = new FastaSequenceIO(outFasta);
		fsio.write(rcs, lineLength);
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-i", "Input fasta", true);
		p.addStringArg("-o", "Output fasta", true);
		p.addIntArg("-l", "Line length", false, 60);
		p.parse(args);
		String in = p.getStringArg("-i");
		String out = p.getStringArg("-o");
		int lineLength = p.getIntArg("-l");
		
		ReverseComplementSequences rcs = new ReverseComplementSequences(in);
		rcs.writeReverseComplements(out, lineLength);
		
		logger.info("All done.");
		
	}

}
