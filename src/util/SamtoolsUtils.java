package util;

import java.util.Collection;

import guttmanlab.core.sequence.FastaFileIOImpl;
import guttmanlab.core.sequence.Sequence;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMSequenceRecord;

public class SamtoolsUtils {
	
	/**
	 * Create a SAM file header for the sequences in a fasta file
	 * @param fastaFile Fasta file
	 * @return SAM header with these sequences and lengths
	 */
	public static final SAMFileHeader createSamHeader(String fastaFile) {
		Collection<Sequence> seqs = new FastaFileIOImpl().readFromFile(fastaFile);
		SAMFileHeader rtrn = new SAMFileHeader();
		for(Sequence seq : seqs) {
			rtrn.addSequence(new SAMSequenceRecord(seq.getName(), seq.getLength()));
		}
		return rtrn;
	}
	
}
