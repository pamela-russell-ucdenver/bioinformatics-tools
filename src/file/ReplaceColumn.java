package file;

import guttmanlab.core.util.CommandLineParser;
import guttmanlab.core.util.StringParser;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;


public class ReplaceColumn {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-i", "Input table", true);
		p.addIntArg("-c", "Column number to replace (0-based)", true);
		p.addStringArg("-r", "File containing replacement entry for each line", true);
		p.addStringArg("-o", "Output table", true);
		p.parse(args);
		String input = p.getStringArg("-i");
		int col = p.getIntArg("-c");
		String newVals = p.getStringArg("-r");
		String output = p.getStringArg("-o");
		
		if(col < 0) {
			throw new IllegalArgumentException("Column number must be >= 0");
		}
		
		FileReader r = new FileReader(input);
		BufferedReader b = new BufferedReader(r);
		FileReader rn = new FileReader(newVals);
		BufferedReader bn = new BufferedReader(rn);
		FileWriter w = new FileWriter(output);
		StringParser s = new StringParser();
		
		while(b.ready()) {
			String line = b.readLine();
			String newVal = bn.readLine();
			s.parse(line);
			if(s.getFieldCount() == 0) {
				w.write("\n");
				continue;
			}
			if(s.getFieldCount() <= col) {
				r.close();
				b.close();
				rn.close();
				bn.close();
				w.close();
				throw new IllegalArgumentException("Line must have at least " + Integer.valueOf(col+1).toString() + " fields:\n" + line);
			}
			String[] fields = s.getStringArray();
			s.parse(newVal);
			if(s.getFieldCount() > 1) {
				r.close();
				b.close();
				rn.close();
				bn.close();
				w.close();
				throw new IllegalArgumentException("Line must contain only one field:\n" + newVal);
			}
			fields[col] = newVal;
			String newLine = fields[0];
			for(int i = 1; i < fields.length; i++) {
				newLine += "\t" + fields[i];
			}
			w.write(newLine + "\n");
		}
		
		if(bn.ready()) {
			r.close();
			b.close();
			rn.close();
			bn.close();
			w.close();
			throw new IllegalArgumentException("Files must have same number of lines");
		}
		
		r.close();
		b.close();
		rn.close();
		bn.close();
		w.close();

	}

}
