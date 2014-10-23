package annotation;

import guttmanlab.core.util.CommandLineParser;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;

import org.apache.log4j.Logger;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Gene;

import broad.pda.annotation.BEDFileParser;

/**
 * @author prussell
 *
 */
public class STARSpliceJunctions {
	
	private static Logger logger = Logger.getLogger(STARSpliceJunctions.class.getName());
	
	private static void writeSpliceJunctions(String inputBed, String outputTable) throws IOException {
		logger.info("Reading genes from file " + inputBed + " and writing splice junctions to " + outputTable + "...");
		Map<String, Collection<Gene>> genesByChr = BEDFileParser.loadDataByChr(new File(inputBed));
		FileWriter w = new FileWriter(outputTable);
		for(String chr : genesByChr.keySet()) {
			for(Gene gene : genesByChr.get(chr)) {
				Collection<? extends Annotation> introns = gene.getIntronSet();
				for(Annotation intron : introns) {
					int oneBasedStart = intron.getStart() + 1;
					int oneBasedEndInclusive = intron.getEnd();
					String strand = intron.getStrand().toString();
					w.write(chr + "\t" + oneBasedStart + "\t" + oneBasedEndInclusive + "\t" + strand + "\n");
				}
			}
		}
		w.close();
		logger.info("Done writing table.");
	}
	
	/**
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-i", "Input bed file", true);
		p.addStringArg("-o", "Output bed file", true);
		p.parse(args);
		String in = p.getStringArg("-i");
		String out = p.getStringArg("-o");
		
		writeSpliceJunctions(in, out);
		
		logger.info("All done.");
		
	}

}
