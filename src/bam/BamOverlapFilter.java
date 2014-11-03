package bam;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;

import org.apache.log4j.Logger;

import net.sf.samtools.BAMFileWriter;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.FeatureCollection;
import guttmanlab.core.util.CommandLineParser;

/**
 * Filter a bam file based on region overlaps, ignoring strand
 * @author prussell
 *
 */
public class BamOverlapFilter {
	
	private String bamFile;
	private static Logger logger = Logger.getLogger(BamOverlapFilter.class.getName());
	
	private BamOverlapFilter(String inputBam) {
		bamFile = inputBam;
	}
	
	private static FeatureCollection<Gene> getExons(FeatureCollection<Gene> genes) {
		FeatureCollection<Gene> rtrn = new FeatureCollection<Gene>(null);
		for(Gene gene : genes) {
			Iterator<SingleInterval> iter = gene.getBlocks();
			while(iter.hasNext()) {
				SingleInterval interval = iter.next();
				rtrn.add(new Gene(interval));
			}
		}
		return rtrn;
	}
	
	private static FeatureCollection<Gene> getSpans(FeatureCollection<Gene> genes) {
		FeatureCollection<Gene> rtrn = new FeatureCollection<Gene>(null);
		for(Gene gene : genes) {
			Gene span = new Gene(new SingleInterval(gene.getReferenceName(), gene.getReferenceStartPosition(), gene.getReferenceEndPosition(), gene.getOrientation()));
			rtrn.add(span);
		}
		return rtrn;
	}
	
	private void writeReadsThatOverlapExons(FeatureCollection<Gene> genes, boolean primaryAlignmentsOnly, String outputBam) {
		writeFilteredFile(getExons(genes), true, primaryAlignmentsOnly, outputBam);
	}
	
	private void writeReadsThatOverlapGeneSpans(FeatureCollection<Gene> genes, boolean primaryAlignmentsOnly, String outputBam) {
		writeFilteredFile(getSpans(genes), true, primaryAlignmentsOnly, outputBam);
	}
	
	private void writeReadsThatDoNotOverlapExons(FeatureCollection<Gene> genes, boolean primaryAlignmentsOnly, String outputBam) {
		writeFilteredFile(getExons(genes), false, primaryAlignmentsOnly, outputBam);
	}
	
	private void writeReadsThatDoNotOverlapGeneSpans(FeatureCollection<Gene> genes, boolean primaryAlignmentsOnly, String outputBam) {
		writeFilteredFile(getSpans(genes), false, primaryAlignmentsOnly, outputBam);
	}
	
	
	private void writeFilteredFile(FeatureCollection<Gene> regions, boolean keepOverlappers, boolean primaryAlignmentsOnly, String outputBam) {
		logger.info("");
		logger.info("Writing to file " + outputBam + "...");
		SAMFileReader reader = new SAMFileReader(new File(bamFile));
		SAMFileHeader header = reader.getFileHeader();
		SAMRecordIterator iter = reader.iterator();
		
		BAMFileWriter writer=new BAMFileWriter(new File(outputBam));
		writer.setHeader(header);
		
		int numRecords = 0;
		int numWritten = 0;
		
		while(iter.hasNext()) {
			SAMRecord record=iter.next();
			numRecords++;
			if(numRecords % 1000000 == 0) {
				logger.info("Wrote " + numWritten + " records out of " + numRecords + ".");
			}
			Annotation mappedRegion = new SingleInterval(record.getReferenceName(), record.getAlignmentStart(), record.getAlignmentEnd(), Strand.BOTH);
			Boolean overlaps = regions.overlaps(mappedRegion);
			if(primaryAlignmentsOnly && record.getNotPrimaryAlignmentFlag()) {
				continue;
			}
			if((overlaps && keepOverlappers) || (!overlaps && !keepOverlappers)) {
				writer.addAlignment(record);
				numWritten++;
			}
		}
		reader.close();
		writer.close();
		logger.info("Done writing file.");
	}
	
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-i", "Input bam", true);
		p.addStringArg("-a", "Bed annotation for overlap filter", true);
		p.addStringArg("-o", "Output filtered bam", true);
		p.addBooleanArg("-ko", "Keep overlappers. If false, remove overlappers.", true);
		p.addBooleanArg("-ex", "Only use exons. If false, use entire gene spans.", true);
		p.addBooleanArg("-pr", "Keep primary alignments only", false, true);
		p.addStringArg("-r", "Reference sequence length file", true);
		p.parse(args);
		String inputBam = p.getStringArg("-i");
		String bed = p.getStringArg("-a");
		String outputBam = p.getStringArg("-o");
		boolean keepOverlappers = p.getBooleanArg("-ko");
		boolean exonsOnly = p.getBooleanArg("-ex");
		String refLengths = p.getStringArg("-r");
		boolean primaryAlignmentsOnly = p.getBooleanArg("-pr");
		
		FeatureCollection<Gene> genes = (FeatureCollection<Gene>) BEDFileIO.loadFromFile(bed, refLengths);
		BamOverlapFilter b = new BamOverlapFilter(inputBam);
		
		if(exonsOnly && keepOverlappers) {
			b.writeReadsThatOverlapExons(genes, primaryAlignmentsOnly, outputBam);
		}
		
		if(exonsOnly && !keepOverlappers) {
			b.writeReadsThatDoNotOverlapExons(genes, primaryAlignmentsOnly, outputBam);
		}
		
		if(!exonsOnly && keepOverlappers) {
			b.writeReadsThatOverlapGeneSpans(genes, primaryAlignmentsOnly, outputBam);
		}
		
		if(!exonsOnly && !keepOverlappers) {
			b.writeReadsThatDoNotOverlapGeneSpans(genes, primaryAlignmentsOnly, outputBam);
		}
		
		logger.info("");
		logger.info("All done.");
		
	}

}
