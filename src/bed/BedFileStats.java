package bed;

import guttmanlab.core.util.CommandLineParser;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;

import org.apache.log4j.Logger;

import nextgen.core.annotation.Gene;

import broad.core.math.Statistics;
import broad.pda.annotation.BEDFileParser;

/**
 * @author prussell
 *
 */
public class BedFileStats {

	private Map<String, Collection<Gene>> genes;
	private static Logger logger = Logger.getLogger(BedFileStats.class.getName());
	
	/**
	 * @param bedFile Bed file of genes
	 * @throws IOException
	 */
	public BedFileStats(String bedFile) throws IOException {
		genes = BEDFileParser.loadDataByChr(new File(bedFile));
	}
	
	/**
	 * Get median size of genes
	 * @return Median gene size
	 */
	public double medianGeneSize() {
		ArrayList<Double> sizes = new ArrayList<Double>();
		for(String chr : genes.keySet()) {
			for(Gene gene : genes.get(chr)) {
				sizes.add(Double.valueOf(gene.getSize()));
			}
		}
		return Statistics.median(sizes);
	}
	
	private void writeBasicStats(String outFile) throws IOException {
		logger.info("Writing basic stats to file " + outFile + "...");
		FileWriter w = new FileWriter(outFile);
		w.write("Median_gene_size\t" + medianGeneSize() + "\n");
		w.close();		
	}
	
	private void writeTableFeatureSize(String outFile) throws IOException {
		logger.info("Writing gene sizes to file " + outFile + "...");
		FileWriter w = new FileWriter(outFile);
		String header = "Gene_name\t";
		header += "Size\t";
		header += "CDS_size\t";
		header += "5UTR_size\t";
		header += "3UTR_size\t";
		w.write(header + "\n");
		for(String chr : genes.keySet()) {
			for(Gene gene : genes.get(chr)) {
				Gene utr5 = gene.get5UTRGene();
				Gene utr3 = gene.get3UTRGene();
				Gene cds = gene.getCDS();
				String line = gene.getName() + "\t";
				line += gene.getSize() + "\t";
				if(cds == null) {
					line += "-\t";
				} else {
					line += cds.getSize() + "\t";
				}
				if(utr5 == null) {
					line += "-\t";
				} else {
					line += utr5.getSize() + "\t";
				}
				if(utr3 == null) {
					line += "-\t";
				} else {
					line += utr3.getSize() + "\t";
				}
				w.write(line + "\n");
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
		p.addStringArg("-g", "Gene bed file", true);
		p.addStringArg("-o", "Output file for basic stats", false, null);
		p.addStringArg("-ot", "Output file for table of gene and feature sizes", false, null);
		p.parse(args);
		String geneFile = p.getStringArg("-g");
		String outFile = p.getStringArg("-o");
		String outTableGeneSize = p.getStringArg("-ot");
		
		BedFileStats b = new BedFileStats(geneFile);
		
		if(outFile != null) {
			b.writeBasicStats(outFile);
		}
		
		if(outTableGeneSize != null) {
			b.writeTableFeatureSize(outTableGeneSize);
		}
		
		logger.info("All done.");
		
	}

}
