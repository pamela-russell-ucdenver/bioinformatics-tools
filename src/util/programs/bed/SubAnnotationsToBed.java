package util.programs.bed;

import guttmanlab.core.util.CommandLineParser;
import guttmanlab.core.util.StringParser;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Map;

import org.apache.log4j.Logger;

import broad.pda.annotation.BEDFileParser;

import nextgen.core.annotation.Gene;

/**
 * @author prussell
 *
 */
public class SubAnnotationsToBed {
	
	private Map<String, Gene> genesByName;
	private static Logger logger = Logger.getLogger(SubAnnotationsToBed.class.getName());
	
	private SubAnnotationsToBed(String geneBed) throws IOException {
		genesByName = BEDFileParser.loadDataByName(new File(geneBed));
	}
	
	private Gene getSubGene(String geneName, int startPosOnTranscript, int endPosOnTranscript) {
		if(!genesByName.containsKey(geneName)) {
			throw new IllegalArgumentException("Gene " + geneName + " not recognized.");
		}
		Gene gene = genesByName.get(geneName);
		int endpoint1 = gene.transcriptToGenomicPosition(startPosOnTranscript);
		int endpoint2 = gene.transcriptToGenomicPosition(endPosOnTranscript);
		int genomicStart = Math.min(endpoint1, endpoint2);
		int genomicEnd = Math.max(endpoint1, endpoint2) + 1;
		return gene.trimAbsolute(genomicStart, genomicEnd);
	}
	
	private void writeSubGeneBed(String coordTable, String outBed) throws IOException {
		logger.info("Writing sub-genes to file " + outBed + "...");
		FileReader r = new FileReader(coordTable);
		BufferedReader b = new BufferedReader(r);
		FileWriter w = new FileWriter(outBed);
		StringParser s = new StringParser();
		while(b.ready()) {
			s.parse(b.readLine());
			String geneName = s.asString(0);
			int start = s.asInt(1);
			int end = s.asInt(2);
			Gene subGene = getSubGene(geneName, start, end);
			String subGeneName = geneName + "_" + start + "_" + end;
			if(s.getFieldCount() == 4) {
				subGeneName += "_" + s.asString(3);
			}
			subGene.setName(subGeneName);
			w.write(subGene.toBED() + "\n");
		}
		r.close();
		b.close();
		w.close();
		logger.info("Done writing file.");
	}
	
	
	/**
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-g", "Bed file of gene annotations", true);
		p.addStringArg("-c", "Table of sub-gene intervals in transcript coordinates. Line format: '<Gene name>   <Transcript start>   <Transcript end>   <Extra tag to add to subgene name (optional)>'", true);
		p.addStringArg("-o", "Output bed file of sub-genes in genomic coordinates", true);
		p.parse(args);
		String geneBed = p.getStringArg("-g");
		String coordTable = p.getStringArg("-c");
		String outBed = p.getStringArg("-o");
		
		SubAnnotationsToBed satb = new SubAnnotationsToBed(geneBed);
		satb.writeSubGeneBed(coordTable, outBed);
		
		logger.info("All done.");
		
	}

}
