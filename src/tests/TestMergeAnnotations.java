package tests;

import general.CommandLineParser;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;

import org.apache.log4j.Logger;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Gene;
import nextgen.core.utils.AnnotationUtils;

import broad.pda.annotation.BEDFileParser;

public class TestMergeAnnotations {
	
	private static Logger logger = Logger.getLogger(TestMergeAnnotations.class.getName());
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-i", "Input bed", true);
		p.addStringArg("-o", "Output bed", true);
		p.parse(args);
		String in = p.getStringArg("-i");
		String out = p.getStringArg("-o");
		
		Map<String, Collection<Gene>> genes = BEDFileParser.loadDataByChr(in);

		FileWriter w = new FileWriter(out);
		
		for(String chr : genes.keySet()) {
			logger.info(chr);
			Collection<Annotation> merged = AnnotationUtils.mergeOverlappingBlocksSingleOrientation(genes.get(chr), false);
			for(Annotation region : merged) {
				w.write(region.toBED() + "\n");
			}
		}
		w.close();
		
		logger.info("All done.");
	}

}
