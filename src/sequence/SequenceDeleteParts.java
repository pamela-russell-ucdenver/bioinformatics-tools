package sequence;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import guttmanlab.core.util.CommandLineParser;

/**
 * Remove specific subsequences from a sequence
 * e.g. remove rRNAs from 45S to get spacers
 * @author prussell
 *
 */
public class SequenceDeleteParts {
	
	private static Sequence removeAll(Sequence origSequence, Collection<Sequence> seqsToRemove) {
		String origSeq = origSequence.getSequenceBases().toUpperCase();
		String currSeq = origSeq;
		for(Sequence seqToRemove : seqsToRemove) {
			String seq = seqToRemove.getSequenceBases().toUpperCase();
			currSeq = currSeq.replaceAll(seq, "");
		}
		Sequence rtrn = new Sequence(origSequence.getId() + "_subsequences_removed");
		rtrn.setSequenceBases(currSeq);
		return rtrn;
	}
	
	private static List<Sequence> removeAll(Collection<Sequence> origSequences, Collection<Sequence> seqsToRemove) {
		List<Sequence> rtrn = new ArrayList<Sequence>();
		for(Sequence origSequence : origSequences) {
			rtrn.add(removeAll(origSequence, seqsToRemove));
		}
		return rtrn;
	}
	
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-i", "Input fasta file of original sequences", true);
		p.addStringArg("-r", "Sequences to delete if they are subsequences of sequences in the input file", true);
		p.addStringArg("-o", "Output fasta file of original sequences with subsequences removed", true);
		p.parse(args);
		String origSeqsFile = p.getStringArg("-i");
		String toDeleteFile = p.getStringArg("-r");
		String outputFile = p.getStringArg("-o");
				
		Collection<Sequence> origSeqs = FastaSequenceIO.loadSequences(new File(origSeqsFile));
		Collection<Sequence> toDelete = FastaSequenceIO.loadSequences(new File(toDeleteFile));
		List<Sequence> modifiedSeqs = removeAll(origSeqs, toDelete);
		
		FastaSequenceIO writer = new FastaSequenceIO(outputFile);
		writer.write(modifiedSeqs);

	}

}
