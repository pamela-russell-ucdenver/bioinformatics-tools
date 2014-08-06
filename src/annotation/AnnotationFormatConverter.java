/**
 * 
 */
package annotation;

import general.CommandLineParser;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;

import nextgen.core.annotation.Gene;

import broad.pda.annotation.BEDFileParser;

/**
 * @author prussell
 *
 */
public class AnnotationFormatConverter {

	/**
	 * Convert bed file to GTF format
	 * @param inputBed Bed file
	 * @param outputGTF Output GTF file
	 * @throws Exception 
	 */
	public static void bedToGTF(String inputBed, String outputGTF) throws Exception {
		Map<String, Collection<Gene>> genes = BEDFileParser.loadDataByChr(new File(inputBed));
		for(String chr : genes.keySet()) {
			int totalSize = 0;
			int totalBlocks = 0;
			for(Gene gene : genes.get(chr)) {
				totalSize += gene.getSize();
				totalBlocks += gene.getNumExons();
			}
			System.err.println("Chr=" + chr + "\ttotal_genes=" + genes.get(chr).size() + "\ttotal_blocks=" + totalBlocks + "\ttotal_size=" + totalSize);
		}
		FileWriter w = new FileWriter(outputGTF);
		for(String chr : genes.keySet()) {
			for(Gene gene : genes.get(chr)) {
				w.write(gene.toGTF("*", gene.getName(), gene.getName()));
			}
		}
		w.close();
	}
	
	/**
	 * Convert GTF file to bed format
	 * @param inputGTF GTF file
	 * @param outputBed Output bed file
	 * @throws IOException
	 */
	public static void gtfToBed(String inputGTF, String outputBed) throws IOException {
		
		// TODO finish
		
	}
	
	
	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-i", "Input file", true);
		p.addStringArg("-o", "Output file", true);
		p.addBooleanArg("-b2g", "Bed format to GTF format (if false, GTF to bed)", true);
		p.parse(args);
		String input = p.getStringArg("-i");
		String output = p.getStringArg("-o");
		boolean b2g = p.getBooleanArg("-b2g");
		
		if(b2g) bedToGTF(input,output);
		//else gtfToBed(input,output);
		
	}

}
