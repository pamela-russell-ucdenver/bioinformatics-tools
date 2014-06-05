package alignment;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;

import broad.core.parser.CommandLineParser;
import broad.core.sequence.Sequence;

import nextgen.core.alignment.SmithWatermanAlignment;
import nextgen.core.utils.FileUtil;
import jaligner.Alignment;

public class SelfDimer {
	
	private static int kmerStart = -1;
	private static int k = -1;
	
	public static Alignment alignReverseComplement(String sequence) {
		SmithWatermanAlignment s = new SmithWatermanAlignment(sequence, SmithWatermanAlignment.DEFAULT_MATCH_SCORE, SmithWatermanAlignment.DEFAULT_MISMATCH_SCORE, Float.MAX_VALUE, Float.MAX_VALUE);
		String rc = Sequence.reverseSequence(sequence);
		return s.align(rc);
	}
	
	public static Alignment alignReverse(String sequence) {
		SmithWatermanAlignment s = new SmithWatermanAlignment(sequence, SmithWatermanAlignment.DEFAULT_MATCH_SCORE, SmithWatermanAlignment.DEFAULT_MISMATCH_SCORE, Float.MAX_VALUE, Float.MAX_VALUE);
		char[] r = new char[sequence.length()];
		for(int i = 0; i < r.length; i++) {
			r[i] = sequence.charAt(sequence.length() - i - 1);
		}
		String reverse = new String(r);
		return s.align(reverse);
	}
	
	public static void printSelfAlignments(Collection<String> sequences, String outFilePrefix) throws IOException {
		FileWriter wt = new FileWriter(outFilePrefix + "_table.out");
		String header = "";
		if(kmerStart >= 0 && k > 0) {
			header += "Kmer\tKmer_reverse\t";
		}
		header += "Sequence\tGC%\tSelf_alignment_length\tSelf_alignment_num_matches\tLongest_perfect_match\tLongest_alignment_single_mismatch";
		wt.write(header + "\n");
		FileWriter wa = new FileWriter(outFilePrefix + "_alignments.out");
		for(String seq : sequences) {
			Alignment align = alignReverseComplement(seq);
			int len = align.getLength();
			int matches = align.getNumberOfMatches();
			int maxContig = align.getLongestPerfectMatch();
			int maxGapped = align.getLongestUngappedAlignment(1);
			Sequence s = new Sequence("");
			s.setSequenceBases(seq);
			float gc = s.gcContent();
			String line = "";
			if(kmerStart >= 0 && k > 0) {
				String kmer = seq.substring(kmerStart, kmerStart + k);
				String kmerRev = Sequence.reverseSequence(kmer);
				line += kmer + "\t" + kmerRev + "\t";
			}
			line += seq + "\t" + gc + "\t" + len + "\t" + matches + "\t" + maxContig + "\t" + maxGapped + "\n";
			wt.write(line);
			wa.write(SmithWatermanAlignment.getFullPrintableAlignment(align) + "\n\n");
		}
		wt.close();
		wa.close();
	}

	public static void printSelfAlignments(String sequenceListFile, String outFilePrefix) throws IOException {
		printSelfAlignments(FileUtil.fileLinesAsList(sequenceListFile), outFilePrefix);
	}
	
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-i", "Input list of sequences", true);
		p.addStringArg("-o", "Output file prefix", true);
		p.addIntArg("-k", "Kmer length within sequence", false, -1);
		p.addIntArg("-ks", "Kmer start within sequence", false, -1);
		p.parse(args);
		String input = p.getStringArg("-i");
		String output = p.getStringArg("-o");
		
		k = p.getIntArg("-k");
		kmerStart = p.getIntArg("-ks");
		
		printSelfAlignments(input, output);
		
	}
	
}
