package fasta;

import guttmanlab.core.util.CommandLineParser;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;

import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;

/**
 * Evenly split a fasta file into smaller files
 * @author prussell
 *
 */
public class FastaSplitter {
	
	private Collection<Sequence> sequences;
	
	private FastaSplitter(String inputFasta) throws IOException {
		sequences = FastaSequenceIO.loadSequences(new File(inputFasta));
	}
	
	private void writeFiles(String outPrefix, int numFiles) throws IOException {
		ArrayList<BufferedWriter> writers = new ArrayList<BufferedWriter>();
		for(int i = 0; i < numFiles; i++) {
			String fileName = outPrefix + "_" + i + ".fa";
			writers.add(new BufferedWriter(new FileWriter(fileName)));
		}
		FastaSequenceIO fsio = new FastaSequenceIO();
		int n = 0;
		for(Sequence seq : sequences) {
			int index = n % numFiles;
			BufferedWriter writer = writers.get(index);
			fsio.write(seq, writer);
			n++;
		}
		for(BufferedWriter writer : writers) {
			writer.close();
		}
	}
	
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-i", "Input fasta", true);
		p.addStringArg("-o", "Output prefix", true);
		p.addIntArg("-n", "Number of files to write", true);
		p.parse(args);
		String input = p.getStringArg("-i");
		String outPrefix = p.getStringArg("-o");
		int numFiles = p.getIntArg("-n");
		
		FastaSplitter fs = new FastaSplitter(input);
		fs.writeFiles(outPrefix, numFiles);
		
		System.out.println("\nAll done.");

	}

}
