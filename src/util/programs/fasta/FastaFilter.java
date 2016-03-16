/**
 * 
 */
package util.programs.fasta;

import guttmanlab.core.util.CommandLineParser;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.TreeSet;

import org.apache.log4j.Logger;

import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;

/**
 * @author prussell
 *
 */
public class FastaFilter {

	private List<? extends Sequence> sequences;
	private static Logger logger = Logger.getLogger(FastaFilter.class.getName());
	
	private FastaFilter(String fastaFile) throws IOException {
		FastaSequenceIO reader = new FastaSequenceIO(fastaFile);
		sequences = reader.loadAll();
	}
	
	private void filterByMinLength(int minLength) {
		logger.info("Removing sequences shorter than " + minLength + " nt.");
		List<Sequence> seqsToRemove = new ArrayList<Sequence>();
		for(Sequence seq : sequences) {
			if(seq.getLength() < minLength) {
				seqsToRemove.add(seq);
			}
		}
		for(Sequence toRemove : seqsToRemove) {
			sequences.remove(toRemove);
		}
	}
	
	private static String getMatch(String longName, Collection<String> shortNames, boolean exactMatch) {
		for(String s : shortNames) {
			if(exactMatch) {
				if(longName.equals(s)) {
					return s;
				}
			} else {
				if(longName.contains(s)) {
					return s;
				}
			}
		}
		return null;
	}
	
	private void keepByName(Collection<String> seqIdsToKeep, boolean exactMatch) {
		logger.info("Keeping sequences by name.");
		List<Sequence> seqsToRemove = new ArrayList<Sequence>();
		Collection<String> notYetFound = new TreeSet<String>();
		notYetFound.addAll(seqIdsToKeep);
		for(Sequence seq : sequences) {
			String id = seq.getId();
			if(getMatch(id, seqIdsToKeep, exactMatch) == null) {
				logger.debug("Removing sequence " + id + " because name is not on list.");
				seqsToRemove.add(seq);
			} else {
				notYetFound.remove(getMatch(id, seqIdsToKeep, exactMatch));
			}
		}
		for(String id : notYetFound) {
			logger.warn("Sequence not found: " + id);
		}
		for(Sequence toRemove : seqsToRemove) {
			sequences.remove(toRemove);
		}

	}
	
	private static Collection<String> linesFromFile(String file) throws IOException {
		FileReader r = new FileReader(file);
		BufferedReader b = new BufferedReader(r);
		Collection<String> rtrn = new ArrayList<String>();
		while(b.ready()) {
			rtrn.add(b.readLine());
		}
		r.close();
		b.close();
		return rtrn;
	}
	
	
	private void write(String outFile) throws IOException {
		logger.info("Writing to " + outFile + ".");
		FastaSequenceIO writer = new FastaSequenceIO(outFile);
		writer.write(sequences);
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {

		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-i", "Input fasta file", true);
		p.addIntArg("--minlen", "Min sequence length", false, 0);
		p.addStringArg("--names", "File of sequence names to keep", false, null);
		p.addStringArg("-o", "Output fasta file", true);
		p.addBooleanArg("-e", "Exact name match for name filter (as opposed to contains)", false, true);
		p.parse(args);
		String inFasta = p.getStringArg("-i");
		String outFasta = p.getStringArg("-o");
		String namesFile = p.getStringArg("--names");
		int minLen = p.getIntArg("--minlen");
		boolean exactMatch = p.getBooleanArg("-e");
		
		FastaFilter ff = new FastaFilter(inFasta);
		if(minLen > 0) {
			ff.filterByMinLength(minLen);
		}
		if(namesFile != null) {
			ff.keepByName(linesFromFile(namesFile), exactMatch);
		}
		
		ff.write(outFasta);
		
		logger.info("");
		logger.info("All done.");
	}

}
