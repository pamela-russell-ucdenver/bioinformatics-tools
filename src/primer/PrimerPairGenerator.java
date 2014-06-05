/**
 * 
 */
package primer;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import broad.core.parser.CommandLineParser;
import broad.core.primer3.PrimerPair;
import broad.core.primer3.PrimerUtils;

/**
 * @author prussell
 *
 */
public class PrimerPairGenerator {

	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		//test6
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-o", "Output file", true);
		p.addIntArg("-n", "Number of primer pairs to generate", true);
		p.addIntArg("-l", "Primer length", true);
		p.addStringArg("-p3", "primer3_core executable", true);
		p.addDoubleArg("-tm", "Optimal Tm", true);
		p.addBooleanArg("-f", "Write out all fields to be read into PrimerPair constructor", false, false);
		p.addStringArg("-p", "File of existing primer pairs to try, as formatted by PrimerPair.getPrimerFieldsAsStringForConstructor()", false, null);
		p.addStringArg("-pi", "Primer ID prefix", false, null);
		p.parse(args);
		String outFile = p.getStringArg("-o");
		int numPrimers = p.getIntArg("-n");
		int primerLength = p.getIntArg("-l");
		String primer3 = p.getStringArg("-p3");
		double optimalMeltingTemp = p.getDoubleArg("-tm");
		boolean full = p.getBooleanArg("-f");
		String primerFile = p.getStringArg("-p");
		BufferedReader primerReader = null;
		if(primerFile != null) {
			FileReader r = new FileReader(primerFile);
			primerReader = new BufferedReader(r);
		}
		String id = p.getStringArg("-pi");
		
		FileWriter w = new FileWriter(outFile);
			if(!full) {
			String header = "left_primer\t";
			header += "right_primer\t";
			header += "left_primer_tm\t";
			header += "right_primer_tm\t";
			header += "primer_pair_penalty\t";
			w.write(header + "\n");
		}
		for(int i=0; i < numPrimers; i++) {
			String thisId = id == null ? null : id + "_" + Integer.valueOf(i+1).toString();
			PrimerPair primer = PrimerUtils.getOneSyntheticPrimerPair(primerLength, primer3, optimalMeltingTemp, primerReader, thisId);
			if(full) {
				w.write(primer.getPrimerFieldsAsStringForConstructor() + "\n");
			} else {
				String line = primer.getLeftPrimer() + "\t";
				line += primer.getRightPrimer() + "\t";
				line += primer.getLeftPrimerTM() + "\t";
				line += primer.getRightPrimerTM() + "\t";
				line += primer.getPrimerPairPenalty() + "\t";
				w.write(line + "\n");
			}
		}
		w.close();
		
		if(primerReader != null) {
			primerReader.close();
		}

	}

}
