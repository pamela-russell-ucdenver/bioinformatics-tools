package rnaseq;

import general.CommandLineParser;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeSet;

import org.apache.log4j.Logger;

import broad.pda.annotation.BEDFileParser;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Gene;
import nextgen.core.coordinatesystem.GenomicSpace;
import nextgen.core.coordinatesystem.TranscriptomeSpace;
import nextgen.core.model.AlignmentModel;
import nextgen.core.model.score.CountScore;

public class RatioIntronExons {
	
	/**
	 * Transcriptome space for mature transcripts
	 * Genomic space for introns
	 * 
	 * For each gene:
	 * Calculate FPKM of mature transcript
	 * Get all introns
	 * Calculate FPKM of each intron (fully contained only)
	 * Write bed file of introns: set score to ratio; change color if ratio is above cutoff
	 */
	
	private Map<String, Collection<Gene>> genes;
	private AlignmentModel matureData;
	private AlignmentModel intronData;
	private double total;
	private static Logger logger = Logger.getLogger(RatioIntronExons.class.getName());
	
	private RatioIntronExons(String bamFile, String geneBed, String chrSizeFile) throws IOException {
		this(bamFile, geneBed, chrSizeFile, -1);
	}
	
	private RatioIntronExons(String bamFile, String geneBed, String chrSizeFile, int totalFragments) throws IOException {
		logger.info("Instantiating...");
		genes = BEDFileParser.loadDataByChr(new File(geneBed));
		TranscriptomeSpace t = new TranscriptomeSpace(genes);
		GenomicSpace g = new GenomicSpace(chrSizeFile);
		matureData = new AlignmentModel(bamFile, t);
		intronData = new AlignmentModel(bamFile, g);
		if(totalFragments > 0) {
			total = totalFragments;
		} else {
			logger.info("Getting total count...");
			total = intronData.getCount();
			logger.info("There are " + total + " total fragments.");
		}
	}
	
	private static double getRPKM(Annotation region, double count, double total) {
		CountScore s = new CountScore(region, count, count, total);
		return s.getRPKM();
	}
	
	private Collection<Annotation> getIntronsWithRatio(Gene gene) {
		double matureCount = matureData.getCount(gene);
		double matureRpkm = getRPKM(gene, matureCount, total);
		Collection<? extends Annotation> introns = gene.getIntronSet();
		Collection<Annotation> rtrn = new TreeSet<Annotation>();
		for(Annotation intron : introns) {
			double intronCount = intronData.getCount(intron, true);
			double intronRpkm = getRPKM(intron, intronCount, total);
			double ratio = intronRpkm / matureRpkm;
			intron.setScore(ratio);
			
			rtrn.add(intron);
		}
		return rtrn;
	}
	
	private void writeIntronBed(String outBed, double maxRatio, boolean keepAboveRatioDifferentColor) throws IOException {
		logger.info("Writing introns with ratio as score to file " + outBed + ".");
		if(keepAboveRatioDifferentColor) {
			logger.info("Introns with ratio >" + maxRatio + " are colored.");
		} else {
			logger.info("Removing introns with ratio > " + maxRatio + ".");
		}
		FileWriter w = new FileWriter(outBed);
		for(String chr : genes.keySet()) {
			logger.info(chr);
			for(Gene gene : genes.get(chr)) {
				Collection<Annotation> introns = getIntronsWithRatio(gene);
				for(Annotation intron : introns) {
					if(intron.getScore() > maxRatio) {
						if(keepAboveRatioDifferentColor) {
							w.write(intron.toBED(255, 0, 0) + "\n");
						}
					}
					else w.write(intron.toBED(0,0,0) + "\n");
				}
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
		p.addStringArg("-b", "Bam file", true);
		p.addStringArg("-g", "Gene bed file", true);
		p.addStringArg("-c", "Chromosome size file", true);
		p.addDoubleArg("-r", "Max intron to mature ratio", true);
		p.addStringArg("-o", "Output bed", true);
		p.addIntArg("-t", "Total number of fragments if known (saves time)", false, -1);
		p.addBooleanArg("-k", "Keep introns with ratio above cutoff and use different color in bed file", false, false);
		p.parse(args);
		String bamFile = p.getStringArg("-b");
		String chrSizeFile = p.getStringArg("-c");
		String bedFile = p.getStringArg("-g");
		double ratio = p.getDoubleArg("-r");
		boolean keep = p.getBooleanArg("-k");
		if(ratio < 0 || ratio > 1) {
			throw new IllegalArgumentException("Ratio must be between 0 and 1.");
		}
		String outBed = p.getStringArg("-o");
		int totalFragments = p.getIntArg("-t");
		
		RatioIntronExons r = new RatioIntronExons(bamFile, bedFile, chrSizeFile, totalFragments);
		r.writeIntronBed(outBed, ratio, keep);
		
		logger.info("");
		logger.info("All done.");
		
	}

}
