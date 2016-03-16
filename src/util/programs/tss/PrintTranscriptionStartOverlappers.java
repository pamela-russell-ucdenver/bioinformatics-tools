package util.programs.tss;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;

import org.apache.log4j.Logger;

import broad.core.parser.CommandLineParser;
import broad.pda.annotation.BEDFileParser;

import net.sf.samtools.util.CloseableIterator;
import nextgen.core.alignment.Alignment;
import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Annotation.Strand;
import nextgen.core.annotation.BasicAnnotation;
import nextgen.core.annotation.Gene;
import nextgen.core.coordinatesystem.GenomicSpace;
import nextgen.core.model.AlignmentModel;

public class PrintTranscriptionStartOverlappers {
	
	private static Logger logger = Logger.getLogger(PrintTranscriptionStartOverlappers.class.getName());
	
	private static Annotation getFivePrimeEnd(Gene gene) {
		String chr = gene.getChr();
		int start = -1;
		if(gene.getOrientation().equals(Strand.POSITIVE)) {
			start = gene.getStart();
		} else if(gene.getOrientation().equals(Strand.NEGATIVE)) {
			start = gene.getEnd() - 1;
		} else {
			throw new IllegalArgumentException("Strand must be known");
		}
		Annotation rtrn = new BasicAnnotation(chr, start, start+1);
		rtrn.setName(gene.getName() + "_five_prime_end");
		return rtrn;
	}
	
	private static void writeOverlappers(AlignmentModel data, Annotation region, FileWriter writer) throws IOException {
		CloseableIterator<Alignment> iter = data.getOverlappingReads(region, false);
		while(iter.hasNext()) {
			Alignment alignment = iter.next();
			writer.write(region.getName() + "\t" + alignment.toSAMRecord().getSAMString());
		}
		iter.close();
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-b", "Bam file", true);
		p.addStringArg("-a", "Bed annotation", true);
		p.addStringArg("-c", "Chromosome size file", true);
		p.addStringArg("-o", "Output file", true);
		p.parse(args);
		String bamFile = p.getStringArg("-b");
		String bedFile = p.getStringArg("-a");
		String chrFile = p.getStringArg("-c");
		String outFile = p.getStringArg("-o");
		
		AlignmentModel data = new AlignmentModel(bamFile, new GenomicSpace(chrFile), false);
		Map<String, Collection<Gene>> genes = BEDFileParser.loadDataByChr(bedFile);
		FileWriter w = new FileWriter(outFile);
		
		logger.info("Writing 5' end overlappers to " + outFile + "...");
		
		for(String chr : genes.keySet()) {
			logger.info(chr);
			for(Gene gene : genes.get(chr)) {
				writeOverlappers(data, getFivePrimeEnd(gene), w);
			}
		}

		w.close();
		
		logger.info("");
		logger.info("All done.");
		
	}

}
