package util.programs.primer;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;

import broad.core.parser.CommandLineParser;
import broad.core.primer3.Primer3Configuration;
import broad.core.primer3.Primer3ConfigurationFactory;
import broad.core.primer3.Primer3IO;
import broad.core.primer3.Primer3SequenceInputTags;
import broad.core.primer3.PrimerPair;
import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;

public class FullLengthPrimerDesigner {

	private static PrimerPair makePrimers(Sequence seq, double minTm, double maxTm, String pathPrimer3core) throws IOException {
		return makePrimers(seq, minTm, maxTm, pathPrimer3core, 0);
	}
	
	private static PrimerPair makePrimers(Sequence seq, double minTm, double maxTm, String pathPrimer3core, int buffer) throws IOException {
		
		Primer3Configuration config = Primer3ConfigurationFactory.getSyntheticConfiguration((maxTm+minTm)/2);				
        config.minGCContent = 20.0;
        config.maxGCContent = 80.0;
        config.optimalPrimerSize = 20;
        config.minPrimerSize = 10;
        config.maxPrimerSize = 35;
		config.minProductSize = Math.max(seq.getLength() - buffer, config.maxPrimerSize);
        config.maxProductSize = seq.getLength();
        config.overLengthPenaltyWeight=0.0001;
        config.underLengthPenaltyWeight=0.0001;
        config.maxMeltingTemp = maxTm;
        config.minMeltingTemp = minTm;
 		
		Primer3IO p3io = new Primer3IO(pathPrimer3core);
		p3io.startPrimer3Communications();

		Primer3SequenceInputTags p3sit = new Primer3SequenceInputTags(seq);	
		p3sit.setPrimerSequenceId(seq.getId());
		Collection<PrimerPair> pp = p3io.findPrimerPair(p3sit, config);
		p3io.endPrimer3Communications();
		
		try {
			return pp.iterator().next();
		} catch(Exception e) {
			return null;
		}
		
	}

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-f", "Fasta file of sequences", true);
		p.addDoubleArg("-mintm", "Min Tm for primers", true);
		p.addDoubleArg("-maxtm", "Max Tm for primers", true);
		p.addStringArg("-o", "Output", true);
		p.addStringArg("-p", "primer3_core executable", true);
		p.addIntArg("-b", "Buffer length, amount that product size can be less than sequence size", false, 0);
		p.parse(args);
		String fasta = p.getStringArg("-f");
		double mintm = p.getDoubleArg("-mintm");
		double maxtm = p.getDoubleArg("-maxtm");
		String output = p.getStringArg("-o");
		String primer3 = p.getStringArg("-p");
		int buffer = p.getIntArg("-b");
		
		Collection<Sequence> seqs = FastaSequenceIO.loadSequences(new File(fasta));
		FileWriter w = new FileWriter(output);
		
		String header = "seq_ID\t";
		header += "left_or_right_primer\t";
		header += "primer_Tm\t";
		header += "product_size\t";
		w.write(header + "\n");
		
		for(Sequence seq : seqs) {
			PrimerPair primer = makePrimers(seq, mintm, maxtm, primer3, buffer);
			if(primer != null) {
				w.write(seq.getId() + "\tLEFT_PRIMER\t" + primer.getLeftPrimer().toUpperCase() + "\t" + primer.getLeftPrimerTM() + "\t" + primer.getProductSize() + "\n");
				w.write(seq.getId() + "\tRIGHT_PRIMER\t" + primer.getRightPrimer().toUpperCase() + "\t" + primer.getRightPrimerTM() + "\t" + primer.getProductSize() + "\n");
			} else {
				w.write(seq.getId() + "\tNO_VALID_PRIMERS\n");
			}
		}
		
		w.close();

	}

}
