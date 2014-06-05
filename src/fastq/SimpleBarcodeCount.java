package fastq;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Map;
import java.util.TreeMap;

import broad.core.parser.CommandLineParser;

/**
 * @author prussell
 *
 */
public class SimpleBarcodeCount {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub

		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-r", "fastq file of reads", true);
		p.addIntArg("-start", "first position of barcode (0-based)", true);
		p.addIntArg("-end", "last position of barcode (0-based)", true);
		p.addStringArg("-o","outfile",true);
		p.parse(args);
		String fastq = p.getStringArg("-r");
		int start = p.getIntArg("-start");
		int end = p.getIntArg("-end");
		String outfile = p.getStringArg("-o");
		
		Map<String,Integer> counts = new TreeMap<String,Integer>();
		
		FileReader r = new FileReader(fastq);
		BufferedReader b = new BufferedReader(r);
		
		int linenum = 0;
		while(b.ready()) {
			
			String line = b.readLine();
			linenum++;
			
			if(linenum % 4 == 2) {
				String barcode = line.substring(start,end+1);
				if(counts.containsKey(barcode)) {
					Integer newcount = Integer.valueOf(counts.get(barcode).intValue() + 1);
					counts.put(barcode, newcount);
				}
				else {
					counts.put(barcode, Integer.valueOf(1));
				}
			}
			
		}
		
		FileWriter w = new FileWriter(outfile);
		for(String barcode : counts.keySet()) {
			w.write(barcode + "\t" + counts.get(barcode) + "\n");
		}
		
		w.close();
		r.close();
		b.close();
		
	}

}
