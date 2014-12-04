package fasta;

import guttmanlab.core.sequence.FastaFileIOImpl;
import guttmanlab.core.sequence.Sequence;
import guttmanlab.core.util.CommandLineParser;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;


public class FastaReplaceNames {
	
	private static int numLines(String file) throws IOException {
		BufferedReader b = new BufferedReader(new FileReader(file));
		int rtrn = 0;
		while(b.ready()) {
			b.readLine();
			rtrn++;
		}
		b.close();
		return rtrn;
	}
	
	private static void replaceNames(String inputFasta, String inputNamesList, String outputFasta) throws IOException {
		FastaFileIOImpl fio = new FastaFileIOImpl();
		Collection<Sequence> seqs = fio.readFromFile(inputFasta);
		if(seqs.size() != numLines(inputNamesList)) {
			throw new IOException("Number of lines in names list does not match number of sequences");
		}
		seqs = null;
		BufferedReader namesReader = new BufferedReader(new FileReader(inputNamesList));
		BufferedReader fastaReader = new BufferedReader(new FileReader(inputFasta));
		FileWriter writer = new FileWriter(outputFasta);
		while(fastaReader.ready()) {
			String line = fastaReader.readLine();
			if(line.startsWith(">")) {
				writer.write(">" + namesReader.readLine() + "\n");
			} else {
				writer.write(line + "\n");
			}
		}
		namesReader.close();
		fastaReader.close();
		writer.close();
	}
	
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-i", "Input fasta", true);
		p.addStringArg("-n", "File containing list of names to use in fasta file. Must have same number of names as sequences in fasta file.", true);
		p.addStringArg("-o", "Output fasta file", true);
		p.parse(args);
		String inputFasta = p.getStringArg("-i");
		String inputNamesList = p.getStringArg("-n");
		String outputFasta = p.getStringArg("-o");
		
		replaceNames(inputFasta, inputNamesList, outputFasta);

	}

}
