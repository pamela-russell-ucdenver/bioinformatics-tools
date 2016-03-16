package util.programs.bed;

import guttmanlab.core.util.CommandLineParser;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import org.apache.log4j.Logger;

import nextgen.core.annotation.Gene;
import nextgen.core.utils.AnnotationUtils;

import broad.pda.annotation.BEDFileParser;



/**
 * @author prussell
 *
 */
public class BedFileCollapseOverlappers {

	private static Logger logger = Logger.getLogger(BedFileCollapseOverlappers.class.getName());
	
	
	/**
	 * @param inputBed
	 * @param outputBed
	 * @throws IOException
	 */
	public static void collapseOverlappersAndWrite(String inputBed, String outputBed, boolean ignoreStrand) throws IOException {
		logger.info("Collapsing overlapping genes in file " + inputBed);
		Map<String, Collection<Gene>> genes = BEDFileParser.loadDataByChr(new File(inputBed));
		Map<String, Collection<Gene>> collapsed = AnnotationUtils.collapseOverlappers(genes, ignoreStrand);
		logger.info("Done collapsing genes. Writing to file " + outputBed);
		FileWriter w = new FileWriter(outputBed);
		for(String chr : collapsed.keySet()) {
			for(Gene gene : collapsed.get(chr)) {
				w.write(gene.toBED() + "\n");
			}
		}
		w.close();
		logger.info("Done writing.");
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-b", "Input bed file", true);
		p.addStringArg("-o", "Output bed file", true);
		p.addBooleanArg("-s", "Ignore strand", false, false);
		p.parse(args);
		String inputBed = p.getStringArg("-b");
		String outputBed = p.getStringArg("-o");
		boolean ignoreStrand = p.getBooleanArg("-s");
		
		collapseOverlappersAndWrite(inputBed, outputBed, ignoreStrand);
		
	}

}
