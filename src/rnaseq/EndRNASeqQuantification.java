package rnaseq;

import general.CommandLineParser;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;

import org.apache.log4j.Logger;

import broad.pda.annotation.BEDFileParser;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Gene;
import nextgen.core.annotation.Annotation.Strand;
import nextgen.core.coordinatesystem.TranscriptomeSpace;
import nextgen.core.model.AlignmentModel;
import nextgen.core.model.score.ScanStatisticScore;

public class EndRNASeqQuantification {
	
	private Map<String, Collection<Gene>> genes;
	AlignmentModel endData5prime;
	AlignmentModel endData3prime;
	private static Logger logger = Logger.getLogger(EndRNASeqQuantification.class.getName());
	
	private EndRNASeqQuantification(String geneBed, String bam5p, String bam3p) throws IOException {
		genes = BEDFileParser.loadDataByChr(new File(geneBed));
		TranscriptomeSpace ts = new TranscriptomeSpace(genes);
		endData5prime = new AlignmentModel(bam5p, ts);
		endData3prime = new AlignmentModel(bam3p, ts);
	}
	
	private Annotation getLeftEnd(Gene gene, int size) {
		if(gene.getSize() < size) {
			return gene;
		}
		int genomicStart = gene.getStart();
		int genomicEnd = gene.transcriptToGenomicPosition(size - 1) + 1;
		return gene.trimAbsolute(genomicStart, genomicEnd);
	}
	
	private Annotation getRightEnd(Gene gene, int size) {
		if(gene.getSize() < size) {
			return gene;
		}
		int genomicStart = gene.transcriptToGenomicPosition(gene.getSize() - size);
		int genomicEnd = gene.getEnd();
		return gene.trimAbsolute(genomicStart, genomicEnd);
	}
	
	private Annotation get5PrimeEnd(Gene gene, int size) {
		if(gene.getOrientation().equals(Strand.UNKNOWN)) {
			throw new IllegalArgumentException("Strand must be known");
		}
		return gene.getOrientation().equals(Strand.POSITIVE) ? getLeftEnd(gene, size) : getRightEnd(gene, size);
	}
	
	private Annotation get3PrimeEnd(Gene gene, int size) {
		if(gene.getOrientation().equals(Strand.UNKNOWN)) {
			throw new IllegalArgumentException("Strand must be known");
		}
		return gene.getOrientation().equals(Strand.POSITIVE) ? getRightEnd(gene, size) : getLeftEnd(gene, size);
	}
	
	private static double getScanPvalue(AlignmentModel model, Gene gene, Annotation window) {
		double geneCount = model.getCount(gene);
		int geneSize = gene.getSize();
		ScanStatisticScore score = new ScanStatisticScore(model, window, geneCount, geneSize, geneCount, geneSize, false);
		return score.getScanPvalue();
	}
	
	private boolean bothEndsSignificant(Gene gene, int size, double pvalCutoff, int countCutoff, int countCutoffOverride) {
		if(countCutoff >= countCutoffOverride) {
			logger.warn("Count threshold to override scan filter should be >= overall count threshold");
		}
		Annotation end5p = get5PrimeEnd(gene, size);
		Annotation end3p = get3PrimeEnd(gene, size);
		double scan5p = getScanPvalue(endData5prime, gene, end5p);
		double scan3p = getScanPvalue(endData3prime, gene, end3p);
		double count5p = endData5prime.getCount(end5p);
		double count3p = endData3prime.getCount(end3p);
		return (scan5p <= pvalCutoff && scan3p <= pvalCutoff && count5p >= countCutoff && count3p >= countCutoff) || (count5p >= countCutoffOverride && count3p >= countCutoffOverride);
	}
	
	private void writeGenesBothEndsSignificant(int size, double pvalCutoff, int countCutoff, int countCutoffOverride, String outBed) throws IOException {
		logger.info("");
		logger.info("Writing significant isoforms to " + outBed + "...");
		FileWriter w = new FileWriter(outBed);
		for(String chr : genes.keySet()) {
			logger.info(chr);
			for(Gene gene : genes.get(chr)) {
				if(bothEndsSignificant(gene, size, pvalCutoff, countCutoff, countCutoffOverride)) {
					w.write(gene.toBED(255, 0, 0) + "\n");
				}
			}
		}
		w.close();
		logger.info("Done writing file.");
	}
	
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-b", "Bed gene annotation", true);
		p.addStringArg("-a5", "Bam file 5' ends", true);
		p.addStringArg("-a3", "Bam file 3' ends", true);
		p.addIntArg("-w", "Size of end window", true);
		p.addDoubleArg("-p", "Scan P-value cutoff", true);
		p.addStringArg("-o", "Output bed for significant isoforms", true);
		p.addIntArg("-c", "Minimum read count in end regions", true);
		p.addIntArg("-co", "Minimum read count in end regions to override scan cutoff", true);
		p.parse(args);
		String geneBed = p.getStringArg("-b");
		String bam5p = p.getStringArg("-a5");
		String bam3p = p.getStringArg("-a3");
		int windowSize = p.getIntArg("-w");
		double pvalCutoff = p.getDoubleArg("-p");
		String outBed = p.getStringArg("-o");
		int countCutoff = p.getIntArg("-c");
		int countCutoffOverride = p.getIntArg("-co");
		
		EndRNASeqQuantification e = new EndRNASeqQuantification(geneBed, bam5p, bam3p);
		e.writeGenesBothEndsSignificant(windowSize, pvalCutoff, countCutoff, countCutoffOverride, outBed);
		
		logger.info("");
		logger.info("All done.");
		
	}
	
}
