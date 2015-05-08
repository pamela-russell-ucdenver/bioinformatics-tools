package barcode;

import guttmanlab.core.util.CommandLineParser;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.TreeSet;

import broad.core.sequence.Sequence;

public class RandomBarcodeDesigner {
		
	private static boolean containsRepeatedBase(String sequence, int repeatLength) {
		char[] a = new char[repeatLength];
		char[] c = new char[repeatLength];
		char[] g = new char[repeatLength];
		char[] t = new char[repeatLength];
		for (int i=0; i<repeatLength; i++) {
			a[i] = 'A';
			c[i] = 'C';
			g[i] = 'G';
			t[i] = 'T';
		}
		String as = new String(a);
		String cs = new String(c);
		String gs = new String(g);
		String ts = new String(t);
		if(sequence.toUpperCase().contains(as)) return true;
		if(sequence.toUpperCase().contains(cs)) return true;
		if(sequence.toUpperCase().contains(gs)) return true;
		if(sequence.toUpperCase().contains(ts)) return true;
		return false;
	}
	
	private static int numMismatches(String s1, String s2) {
		if(s1.length() != s2.length()) {
			throw new IllegalArgumentException("Strings must have same length");
		}
		char[] s1c = s1.toUpperCase().toCharArray();
		char[] s2c = s2.toUpperCase().toCharArray();
		int rtrn = 0;
		for(int i = 0; i < s1c.length; i++) {
			if(s1c[i] != s2c[i]) rtrn++;
		}
		return rtrn;
	}
	
	private static int minPairwiseMismatches(Collection<String> sequences) {
		int rtrn = Integer.MAX_VALUE;
		for(String seq1 : sequences) {
			for(String seq2 : sequences) {
				if(seq1.equals(seq2)) {
					continue;
				}
				int mm = numMismatches(seq1, seq2);
				if(mm < rtrn) {
					rtrn = mm;
				}
			}
		}
		return rtrn;
	}
	
	private static boolean matchesExisting(String sequence, Collection<String> existingSeqs, int maxMismatches) {
		for(String existing : existingSeqs) {
			if(numMismatches(sequence, existing) <= maxMismatches) {
				return true;
			}
		}
		return false;
	}
	
	/**
	 * Generate a seq of sequences that are pairwise different by a minimum number of mismatches
	 * Also no sequence can have a run of the same base longer than a max length
	 * @param numToGenerate Number of sequences to generate
	 * @param seqLength Length of sequences to generate
	 * @param minMismatches Min required number of pairwise mismatches between all sequences
	 * @param maxRepeatLength Max allowable length of a run of a single base within sequences
	 * @return The collection of generated sequences
	 */
	private static Collection<String> generateNonMatchingSeqsNoRepeatedBases(int numToGenerate, int seqLength, int minMismatches, int maxRepeatLength) {
		int numFailed = 0;
		Collection<String> rtrn = new TreeSet<String>();
		while(rtrn.size() < numToGenerate) {
			if(numFailed > 10000000) {
				break;
			}
			String randSeq = Sequence.generateRandomSequence(seqLength);
			if(containsRepeatedBase(randSeq, maxRepeatLength + 1)) {
				numFailed++;
				continue;
			}
			if(matchesExisting(randSeq, rtrn, minMismatches - 1)) {
				numFailed++;
				continue;
			}
			rtrn.add(randSeq);
		}
		return rtrn;
	}
	
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addIntArg("-n", "Number of barcodes to generate", true);
		p.addIntArg("-l", "Barcode length", true);
		p.addIntArg("-m", "Min pairwise mismatches", true);
		p.addIntArg("-r", "Max length of run of single nucleotide", true);
		p.addStringArg("-o", "Output file", true);
		p.addBooleanArg("-rc", "Also print reverse complement", true);
		p.parse(args);
		int num = p.getIntArg("-n");
		int len = p.getIntArg("-l");
		int mismatch = p.getIntArg("-m");
		int rep = p.getIntArg("-r");
		String out = p.getStringArg("-o");
		boolean revcomp = p.getBooleanArg("-rc");
		
		Collection<String> barcodes = generateNonMatchingSeqsNoRepeatedBases(num, len, mismatch, rep);
		
		FileWriter w = new FileWriter(out);
		for(String barcode : barcodes) {
			String line = barcode;
			if(revcomp) line += "\t" + Sequence.reverseSequence(barcode);
			w.write(line + "\n");
		}		
		w.close();
	}
	
}
