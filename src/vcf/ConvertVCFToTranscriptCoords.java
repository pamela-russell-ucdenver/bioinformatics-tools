package vcf;

import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotationcollection.FeatureCollection;
import guttmanlab.core.util.CommandLineParser;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;

import org.apache.log4j.Logger;

import net.sf.samtools.util.CloseableIterator;
import annotation.OverlapUtils;

/**
 * Convert the coordinates in a VCF file to transcriptome space and reverse complement the allele sequences where applicable
 * There can be multiple transcript conversions for a single VCF record in genome coordinates
 * @author prussell
 *
 */
public class ConvertVCFToTranscriptCoords {
	
	private OverlapUtils overlapUtils;
	private static Logger logger = Logger.getLogger(ConvertVCFToTranscriptCoords.class.getName());
	
	/**
	 * @param geneBedFile Gene annotation in bed format
	 * @param referenceSizes Reference genome size file. Line format: chr   size
	 * @throws IOException
	 */
	private ConvertVCFToTranscriptCoords(String geneBedFile, String referenceSizes) throws IOException {
		overlapUtils = new OverlapUtils(geneBedFile, referenceSizes);
	}
	
	/**
	 * Get the set of converted records for one VCF record
	 * @param record The record
	 * @return Converted records for each transcript overlapping the position
	 */
	private Collection<VCFRecord> getConvertedRecords(VCFRecord record) {
		String chr = record.getChrom();
		int pos = record.getZeroBasedPos();
		FeatureCollection<Gene> geneOverlappers = overlapUtils.getOverlappers(chr, pos);
		CloseableIterator<Gene> iter = geneOverlappers.sortedIterator();
		Collection<VCFRecord> rtrn = new ArrayList<VCFRecord>();
		while(iter.hasNext()) {
			rtrn.add(VCFRecord.convertToTranscriptCoords(record, iter.next()));
		}
		iter.close();
		return rtrn;
	}
	
	/**
	 * Write entire converted VCF file
	 * @param inputVCF Input VCF file
	 * @param outputVCF File to write
	 * @throws IOException
	 */
	private void writeConvertedFile(String inputVCF, String outputVCF) throws IOException {
		logger.info("");
		logger.info("Converting records from " + inputVCF + " to transcript coordinates and writing to " + outputVCF + "...");
		BufferedReader reader = new BufferedReader(new FileReader(inputVCF));
		FileWriter writer = new FileWriter(outputVCF);
		int numDone = 0;
		while(reader.ready()) {
			numDone++;
			if(numDone % 100000 == 0) {
				logger.info("Finished " + numDone + " records.");
			}
			String line = reader.readLine();
			if(line.startsWith("#")) {
				continue;
			}
			VCFRecord record = new VCFRecord(line);
			Collection<VCFRecord> convertedRecords = getConvertedRecords(record);
			for(VCFRecord converted : convertedRecords) {
				writer.write(converted.toString() + "\n");
			}
		}
		reader.close();
		writer.close();
		logger.info("Done writing converted file.");
	}
	
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-g", "Gene annotation bed file", true);
		p.addStringArg("-c", "Chromosome size file", true);
		p.addStringArg("-i", "Input VCF file", true);
		p.addStringArg("-o", "Output VCF file", true);
		p.parse(args);
		String geneBed = p.getStringArg("-g");
		String chrFile = p.getStringArg("-c");
		String inputVCF = p.getStringArg("-i");
		String outputVCF = p.getStringArg("-o");
		
		ConvertVCFToTranscriptCoords c = new ConvertVCFToTranscriptCoords(geneBed, chrFile);
		c.writeConvertedFile(inputVCF, outputVCF);
		
		logger.info("");
		logger.info("All done.");
		
	}
	
}
