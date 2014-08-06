/**
 * 
 */
package motif;

import general.CommandLineParser;
import general.StringParser;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

import org.apache.log4j.Logger;

import broad.pda.annotation.BEDFileParser;

import nextgen.core.annotation.Gene;

/**
 * @author prussell
 *
 */
public class Fimo2Bed {

	protected Map<String, Gene> genes;
	private Map<String, Collection<MotifOccurrence>> motifs;
	protected static Logger logger = Logger.getLogger(Fimo2Bed.class.getName());
	
	private Fimo2Bed(String fimoTxtFile, String bedFile) throws IOException {
		Map<String, Collection<Gene>> genesByChr = BEDFileParser.loadDataByChr(new File(bedFile));
		genes = new TreeMap<String, Gene>();
		for(String chr : genesByChr.keySet()) {
			for(Gene gene : genesByChr.get(chr)) {
				genes.put(gene.getName(), gene);
			}
		}
		motifs = new TreeMap<String, Collection<MotifOccurrence>>();
		FileReader r = new FileReader(fimoTxtFile);
		BufferedReader b = new BufferedReader(r);
		StringParser s = new StringParser();
		while(b.ready()) {
			String line = b.readLine();
			if(line.contains("#")) {
				continue;
			}
			s.parse(line);
			String motif = s.asString(0);
			if(!motifs.containsKey(motif)) {
				Collection<MotifOccurrence> motifOccurrences = new ArrayList<MotifOccurrence>();
				motifs.put(motif, motifOccurrences);
			}
			String geneName = s.asString(1);
			int geneStart = s.asInt(2);
			int geneEnd = s.asInt(3);
			@SuppressWarnings("unused")
			double pVal = s.asDouble(6);
			double qVal = s.asDouble(7);
			String sequence = s.asString(8);
			try {
				MotifOccurrence motifOccurrence = new MotifOccurrence(geneName, geneStart, geneEnd, qVal, sequence);
				motifs.get(motif).add(motifOccurrence);
			} catch (Exception e) {
				logger.warn("Couldn't find gene " + geneName);
				continue;
			}
		}
		r.close();
		b.close();
	}
	
	private void writeBedTracks(String outputDir, String identifier) throws IOException {
		logger.info("Writing bed tracks in directory " + outputDir);
		File outDir = new File(outputDir);
		@SuppressWarnings("unused")
		boolean madeDir = outDir.mkdir();
		if(!outDir.exists()) {
			throw new IOException("Could not create directory " + outputDir);
		}
		for(String motif : motifs.keySet()) {
			logger.info(motif);
			String bedFile = outputDir + "/" + identifier + "_" + motif + ".bed";
			FileWriter w = new FileWriter(bedFile);
			for(MotifOccurrence motifOccurrence : motifs.get(motif)) {
				Gene motifGenomicCoords = motifOccurrence.getMotifGenomicCoords();
				motifGenomicCoords.setName(motifOccurrence.toString());
				w.write(motifGenomicCoords.toBED(53, 146, 192) + "\n");
			}
			w.close();
		}
	}
	
	private class MotifOccurrence {

		public MotifOccurrence(String geneName, int start, int end, double qVal, String seq) {
			if(!genes.containsKey(geneName)) {
				throw new IllegalArgumentException("Gene " + geneName + " not recognized.");
			}
			gene = genes.get(geneName);
			geneStart = start;
			geneEnd = end;
			qValue = qVal;
			sequence = seq;
		}
		
		public Gene getMotifGenomicCoords() {
			/*if(gene.getOrientation().equals(Strand.POSITIVE)) {
				return gene.trimAbsolute(gene.getReferenceCoordinateAtPosition(geneStart-1), gene.getReferenceCoordinateAtPosition(geneEnd));
			}
			if(gene.getOrientation().equals(Strand.NEGATIVE)) {
				return gene.trimAbsolute(gene.getReferenceCoordinateAtPosition(geneEnd), gene.getReferenceCoordinateAtPosition(geneStart-1));
			}
			throw new IllegalStateException("Gene must have orientation (" + gene.getName() + ").");*/
			int deltaStart = geneStart - 1;
			int deltaEnd = gene.getSize() - geneEnd;
			Gene copy = gene.copy();
			logger.info("Gene size " + gene.getSize() + " gene start " + geneStart + " gene end " + geneEnd + " delta start " + deltaStart + " delta end " + deltaEnd);
			copy.trim(deltaStart, deltaEnd);
			return copy;
		}
		
		@Override
		public String toString() {
			return sequence + "_" + qValue;
		}
		
		private Gene gene;
		private String sequence;
		private int geneStart;
		private int geneEnd;
		private double qValue;
		
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-b", "Bed annotation", true);
		p.addStringArg("-f", "fimo.txt file", true);
		p.addStringArg("-id", "Identifier to prepend to bed file names", false, "");
		p.addStringArg("-o", "Output directory", true);
		p.parse(args);
		String bedFile = p.getStringArg("-b");
		String fimoFile = p.getStringArg("-f");
		String outDir = p.getStringArg("-o");
		String id = p.getStringArg("-id");
		
		Fimo2Bed ftb = new Fimo2Bed(fimoFile, bedFile);
		ftb.writeBedTracks(outDir, id);
		

	}

}
