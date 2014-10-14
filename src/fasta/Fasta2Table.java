package fasta;

import general.CommandLineParser;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;

public class Fasta2Table {
	
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-i", "Input fasta", true);
		p.addStringArg("-o", "Output table", true);
		p.parse(args);
		String input = p.getStringArg("-i");
		String output = p.getStringArg("-o");
		
		FastaSequenceIO fsio = new FastaSequenceIO(new File(input));
		List<Sequence> seqs = fsio.loadAll();
		FileWriter w = new FileWriter(output);
		for(Sequence seq : seqs) {
			w.write(seq.getId() + "\t" + seq.getSequenceBases() + "\n");
		}
		w.close();
		
	}
	
}
