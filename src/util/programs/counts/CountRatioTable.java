package util.programs.counts;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Map;
import java.util.TreeMap;

import org.apache.log4j.Logger;

import net.sf.samtools.util.CloseableIterator;
import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.MappedFragment;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.AnnotationCollection;
import guttmanlab.core.annotationcollection.BAMFragmentCollectionFactory;
import guttmanlab.core.util.CommandLineParser;
import guttmanlab.core.util.CountLogger;

public class CountRatioTable {
	
	private String controlName;
	private Map<String, AnnotationCollection<? extends MappedFragment>> data;
	private Map<String, Double> totalCounts;
	private static Logger logger = Logger.getLogger(CountRatioTable.class.getName());
	
	private CountRatioTable(String controlBam, String bamFileList, boolean forceSingleEnd) throws IOException {
		
		logger.info("");
		logger.info("Instantiating...");
		
		controlName = makeSampleName(controlBam);
		
		data = new TreeMap<String, AnnotationCollection<? extends MappedFragment>>();
		BufferedReader reader = new BufferedReader(new FileReader(bamFileList));
		while(reader.ready()) {
			String bam = reader.readLine();
			String sampleName = makeSampleName(bam);
			data.put(sampleName, BAMFragmentCollectionFactory.createFromBam(bam, forceSingleEnd));
		}
		reader.close();
		
		totalCounts = new TreeMap<String, Double>();
		logger.info("Calculating total reads in sample " + controlName);
		totalCounts.put(controlName, Double.valueOf(data.get(controlName).getNumAnnotations()));
		logger.info("TOTAL_COUNT\t" + controlName + "\t" + totalCounts.get(controlName));
		for(String sample : data.keySet()) {
			logger.info("Calculating total reads in sample " + sample);
			totalCounts.put(sample, Double.valueOf(data.get(sample).getNumAnnotations()));
			logger.info("TOTAL_COUNT\t" + sample + "\t" + totalCounts.get(sample));
		}
	
		logger.info("Done instantiating.");
		logger.info("");
		
	}
	
	private static String makeSampleName(String bamFileName) {
		return bamFileName.replaceAll(".bam", "");
	}
	
	private double getCount(Annotation region, String sampleName) {
		return data.get(sampleName).numOverlappers(region, false);
	}
	
	private double getRPKM(Annotation region, String sampleName) {
		AnnotationCollection<? extends MappedFragment> sampleData = data.get(sampleName);
		double count = sampleData.numOverlappers(region, false);
		double regionSize = region.size();
		return (1000000000.0) * count / (totalCounts.get(sampleName) * regionSize);
	}
	
	private double getRpkmRatioVsControl(Annotation region, String sampleName) {
		return getRPKM(region, sampleName) / getRPKM(region, controlName);
	}
	
	private String getHeaderCount() {
		String rtrn = "";
		for(String sample : data.keySet()) {
			rtrn += "count_" + sample + "\t";
		}
		return rtrn;
	}
	
	private String getLineCount(Annotation region) {
		String rtrn = "";
		for(String sample : data.keySet()) {
			rtrn += getCount(region, sample) + "\t";
		}
		return rtrn;
	}
	
	private String getHeaderRPKM() {
		String rtrn = "";
		for(String sample : data.keySet()) {
			rtrn += "rpkm_" + sample + "\t";
		}
		return rtrn;
	}
	
	private String getLineRPKM(Annotation region) {
		String rtrn = "";
		for(String sample : data.keySet()) {
			rtrn += getRPKM(region, sample) + "\t";
		}
		return rtrn;
	}
	
	private String getHeaderRpkmRatioVsControl() {
		String rtrn = "";
		for(String sample : data.keySet()) {
			rtrn += "rpkm_ratio_" + sample + "_" + controlName + "\t";
		}
		return rtrn;
	}
	
	private String getLineRpkmRatioVsControl(Annotation region) {
		String rtrn = "";
		for(String sample : data.keySet()) {
			rtrn += getRpkmRatioVsControl(region, sample) + "\t";
		}
		return rtrn;
	}
	
	private void writeTable(String bedFile, String referenceSizes, String outFile) throws IOException {
		
		logger.info("Writing table to " + outFile + "...");
		
		AnnotationCollection<Gene> genes = BEDFileIO.loadFromFile(bedFile, referenceSizes);
		
		FileWriter fw = new FileWriter(outFile);
		
		String header = "gene\t";
		header += "gene_length\t";
		header += getHeaderRpkmRatioVsControl();
		
		fw.write(header + "\n");
		
		CountLogger cl = new CountLogger(genes.getNumAnnotations(), 10);
		
		CloseableIterator<Gene> geneIter = genes.sortedIterator();
		while(geneIter.hasNext()) {
			cl.advance();
			Gene gene = geneIter.next();
			int length = gene.size();
			String line = gene.getName() + "\t";
			line += length + "\t";
			line += getLineRpkmRatioVsControl(gene);
			fw.write(line + "\n");
		}
		
		fw.close();
		
		logger.info("");
		logger.info("Done writing table.");
		
	}
	
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-c", "Control bam file", true);
		p.addStringArg("-b", "Bam file list", true);
		p.addStringArg("-g", "Gene bed file", true);
		p.addStringArg("-o", "Output table", true);
		p.addStringArg("-r", "Reference size file", true);
		p.addBooleanArg("-s", "Force single end", true);
		p.parse(args);
		String controlBam = p.getStringArg("-c");
		String referenceSizes = p.getStringArg("-r");
		String bamList = p.getStringArg("-b");
		String outTable = p.getStringArg("-o");
		boolean forceSingleEnd = p.getBooleanArg("-s");
		String geneBed = p.getStringArg("-g");
		
		CountRatioTable crt = new CountRatioTable(controlBam, bamList, forceSingleEnd);
		crt.writeTable(geneBed, referenceSizes, outTable);
		
		logger.info("");
		logger.info("All done.");
		
	}

}
