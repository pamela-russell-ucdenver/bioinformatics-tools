/**
 * 
 */
package annotation;

import guttmanlab.core.util.CommandLineParser;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;

import org.apache.log4j.Logger;

import nextgen.core.annotation.Gene;
import nextgen.core.feature.GeneWindow;

import broad.pda.annotation.BEDFileParser;

/**
 * @author prussell
 *
 */
public class WindowWriter {

	private Collection<Gene> genes;
	private static Logger logger = Logger.getLogger(WindowWriter.class.getName());
	
	/**
	 * @param genesBed
	 * @throws IOException
	 */
	public WindowWriter(String genesBed) throws IOException {
		genes = BEDFileParser.loadData(new File(genesBed));
	}
	
	/**
	 * @param windowSize
	 * @param stepSize
	 * @param outFile
	 * @throws IOException
	 */
	public void writeWindowsToFile(int windowSize, int stepSize, String outFile) throws IOException {
		logger.info("Writing windows to file " + outFile + ". Window size = " + windowSize + " . Step size = " + stepSize + ".");
		FileWriter w = new FileWriter(outFile);
		for(Gene gene : genes) {
			Collection<GeneWindow> windows = gene.getWindows(windowSize, stepSize, 0);
			for(GeneWindow window : windows) {
				window.setName(gene.getName() + ":" + window.getChr() + ":" + window.getStart() + "-" + window.getEnd());
				w.write(window.toBED() + "\n");
			}
		}
		w.close();
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-b", "Bed file of genes", true);
		p.addIntArg("-w", "Window size", true);
		p.addIntArg("-s", "Step size", true);
		p.addStringArg("-o", "Output bed file of windows", true);
		p.parse(args);
		String bedFile = p.getStringArg("-b");
		int windowSize = p.getIntArg("-w");
		int stepSize = p.getIntArg("-s");
		String outFile = p.getStringArg("-o");
		
		WindowWriter w = new WindowWriter(bedFile);
		w.writeWindowsToFile(windowSize, stepSize, outFile);

	}

}
