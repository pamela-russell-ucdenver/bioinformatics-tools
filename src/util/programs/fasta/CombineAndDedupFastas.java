package util.programs.fasta;

import guttmanlab.core.sequence.FastaFileIOImpl;
import guttmanlab.core.sequence.Sequence;
import guttmanlab.core.util.CommandLineParser;

import java.util.Collection;
import java.util.HashSet;

import org.apache.log4j.Logger;

public class CombineAndDedupFastas {
	
	private static Logger logger = Logger.getLogger(CombineAndDedupFastas.class.getName());
	
	private static Collection<Sequence> combine(String fasta1, String fasta2, String fasta3) {
		Collection<Sequence> rtrn = new HashSet<Sequence>();
		Collection<Sequence> col1 = new FastaFileIOImpl().readFromFile(fasta1);
		Collection<Sequence> col2 = new FastaFileIOImpl().readFromFile(fasta2);
		Collection<Sequence> col1upper = new HashSet<Sequence>();
		for(Sequence seq : col1) {
			Sequence seqUpper = new Sequence(seq.getName(), seq.getSequenceBases().toUpperCase());
			col1upper.add(seqUpper);
		}
		rtrn.addAll(col1upper);
		Collection<Sequence> col2upper = new HashSet<Sequence>();
		for(Sequence seq : col2) {
			Sequence seqUpper = new Sequence(seq.getName(), seq.getSequenceBases().toUpperCase());
			col2upper.add(seqUpper);
		}
		rtrn.addAll(col2upper);
		if(fasta3 != null) {
			Collection<Sequence> col3 = new FastaFileIOImpl().readFromFile(fasta3);
			Collection<Sequence> col3upper = new HashSet<Sequence>();
			for(Sequence seq : col3) {
				Sequence seqUpper = new Sequence(seq.getName(), seq.getSequenceBases().toUpperCase());
				col3upper.add(seqUpper);
			}
			rtrn.addAll(col3upper);
		}
		logger.info("Returning " + rtrn.size() + " combined sequences.");
		return rtrn;
	}
	
	
	
	private static void combineAndWrite(String fasta1, String fasta2, String fasta3, String outFasta) {
		Collection<Sequence> combined = combine(fasta1, fasta2, fasta3);
		new FastaFileIOImpl().writeToFile(combined, outFasta, 100);
	}
	
	public static void main(String[] args) {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-f1", "Fasta file 1", true);
		p.addStringArg("-f2", "Fasta file 2", true);
		p.addStringArg("-f3", "Fasta file 3", false, null);
		p.addStringArg("-o", "Output fasta", true);
		p.parse(args);
		String fasta1 = p.getStringArg("-f1");
		String fasta2 = p.getStringArg("-f2");
		String fasta3 = p.getStringArg("-f3");
		String outFasta = p.getStringArg("-o");
		
		combineAndWrite(fasta1, fasta2, fasta3, outFasta);
	}

}
