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
import java.util.TreeSet;

import nextgen.core.annotation.Gene;


import broad.pda.annotation.BEDFileParser;

/**
 * @author prussell
 *
 */
public class GeneCoordinateConverter {

	private Map<String, Collection<Gene>> genesByChr;
	private Map<String, Gene> genesByName;
	
	/**
	 * @param geneBedFile Bed file of genes in genomic coordinates
	 * @throws IOException
	 */
	public GeneCoordinateConverter(String geneBedFile) throws IOException {
		this.genesByChr = BEDFileParser.loadDataByChr(new File(geneBedFile));
		this.genesByName = BEDFileParser.loadDataByName(new File(geneBedFile));
	}
	
	/**
	 * Map a set of intervals from genomic coordinates to transcriptome space
	 * @param bedIntervals The intervals in genomic coordinates
	 * @param outIntervalsTranscriptCoords The output bed file
	 * @param endMinus1 Subtract 1 from end coordinate of interval
	 * @throws IOException
	 */
	public void genomicIntervalsToTranscriptCoordinates(String bedIntervals, String outIntervalsTranscriptCoords, boolean endMinus1) throws IOException {
		Map<String, Collection<Gene>> intervals = BEDFileParser.loadDataByChr(new File(bedIntervals));
		FileWriter w = new FileWriter(outIntervalsTranscriptCoords);
		for(String chr : this.genesByChr.keySet()) {
			if(!intervals.containsKey(chr)) continue;
			for(Gene interval : intervals.get(chr)) {
				for(Gene gene : this.genesByChr.get(chr)) {
					if(gene.overlapsExon(interval)) {
						int geneCoordStart = gene.genomicToTranscriptPosition(interval.getStart());
						int end = interval.getEnd();
						if(endMinus1) end -= 1;
						int geneCoordEnd = gene.genomicToTranscriptPosition(end);
						if(geneCoordStart < 0 || geneCoordEnd < 0) continue;
						w.write(gene.getName() + "\t" + Math.min(geneCoordStart, geneCoordEnd) + "\t" + Math.max(geneCoordStart, geneCoordEnd) + "\n");
					}
				}
			}
		}
		w.close();
	}
	
	/**
	 * Map a set of intervals to transcriptome space in the same gene only
	 * @param bedIntervals The intervals in genomic coordinates
	 * @param outIntervalsTranscriptCoords Output bed file
	 * @param endMinus1 Subtract 1 from end coordinate of interval
	 * @throws IOException
	 */
	public void genomicSubtranscriptIntervalsToTranscriptCoordinates(String bedIntervals, String outIntervalsTranscriptCoords, boolean endMinus1) throws IOException {
		Map<String,Gene> intervals = BEDFileParser.loadDataByName(new File(bedIntervals));
		FileWriter w = new FileWriter(outIntervalsTranscriptCoords);
		for(String name : intervals.keySet()) {
			Gene interval = intervals.get(name);
			if(!this.genesByName.containsKey(name)) {
				System.err.println("Gene " + name + " not found.");
				continue;
			}
			Gene gene = this.genesByName.get(name);
			if(gene.overlaps(interval)) {
				int geneCoordStart = gene.genomicToTranscriptPosition(interval.getStart());
				int end = interval.getEnd();
				if(endMinus1) end -= 1;
				int geneCoordEnd = gene.genomicToTranscriptPosition(end);
				if(geneCoordStart < 0 || geneCoordEnd < 0) {
					System.err.println("For interval " + name + " got start " + geneCoordStart + " and end " + geneCoordEnd);
					continue;
				}
				w.write(gene.getName() + "\t" + Math.min(geneCoordStart, geneCoordEnd) + "\t" + Math.max(geneCoordStart, geneCoordEnd) + "\n");
			} else {
				System.err.println("Interval " + name + " does not overlap gene " + gene.getName());
			}
		}	
		w.close();
	}

	
	/**
	 * Map a set of intervals in transcriptome space to genome space
	 * @param bedIntervals The intervals in transcriptome space, with the gene name as chromosome
	 * @param outGenomicCoords The output bed file
	 * @throws IOException
	 */
	public void transcriptIntervalsToGenomicCoordinates(String bedIntervals, String outGenomicCoords) throws IOException {
		Map<String, Collection<Gene>> intervals = BEDFileParser.loadDataByChr(new File(bedIntervals));
		FileWriter w = new FileWriter(outGenomicCoords);
		TreeSet<String> genesNotFound = new TreeSet<String>();
		for(String chr : intervals.keySet()) {
			for(Gene interval : intervals.get(chr)) {
				String geneName = interval.getChr();
				if(!this.genesByName.containsKey(geneName)) {
					if(!genesNotFound.contains(geneName)) System.err.println("Warning: gene " + geneName + " not found.");
					genesNotFound.add(geneName);
					continue;
				}
				Gene theGene = this.genesByName.get(geneName);
				int geneSize = theGene.getSize();
				int intervalStart = interval.getStart();
				int intervalEnd = interval.getEnd();
				if(intervalStart >= geneSize) {
					System.err.println("WARNING interval start is after gene end. Start=" + intervalStart + " gene size=" + geneSize);
					continue;
				}
				if(intervalEnd >= geneSize) {
					System.err.println("WARNING interval end is after gene end. Trimming interval back to size of gene. End=" + intervalEnd + " gene size =" + geneSize);
					intervalEnd = geneSize-1;
				}
				Gene overlap = theGene.transcriptToGenomicPosition(intervalStart, intervalEnd);
				overlap.setName(interval.getName());
				overlap.setBedScore(interval.getBedScore());
				w.write(overlap.toBED() + "\n");
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
		p.addStringArg("-f", "Interval bed file", true);
		p.addStringArg("-o", "Output bed file", true);
		p.addBooleanArg("-t2g", "Convert intervals from transcriptome space to genome space", true);
		p.parse(args);
		String geneBed = p.getStringArg("-g");
		String featureBed = p.getStringArg("-f");
		String out = p.getStringArg("-o");
		boolean t2g = p.getBooleanArg("-t2g");
		
		GeneCoordinateConverter gcc = new GeneCoordinateConverter(geneBed);
		if(t2g) gcc.transcriptIntervalsToGenomicCoordinates(featureBed, out);
		else gcc.genomicSubtranscriptIntervalsToTranscriptCoordinates(featureBed, out, true);
		
	}

}
