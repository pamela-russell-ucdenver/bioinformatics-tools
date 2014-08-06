package rnaseq;

import general.CommandLineParser;

import java.io.File;
import java.io.IOException;

import nextgen.core.annotation.Gene;
import nextgen.core.coordinatesystem.GenomicSpace;
import nextgen.core.coordinatesystem.TranscriptomeSpace;
import nextgen.core.model.AlignmentModel;
import broad.pda.annotation.BEDFileParser;

/**
 * @author prussell
 *
 */
public class TestReadCounting {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {

		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-b", "Bam file", true);
		p.addStringArg("-g", "Bed gene annotation", true);
		p.addStringArg("-c", "Chromosome size file", true);
		p.addStringArg("-chr", "Chromosome name", true);
		p.addIntArg("-s", "Start position", true);
		p.addIntArg("-e", "End position", true);
		p.addBooleanArg("-pe", "Use fragments", true);
		p.parse(args);
		String bamFile = p.getStringArg("-b");
		String bedFile = p.getStringArg("-g");
		String chrSizeFile = p.getStringArg("-c");
		String chrName = p.getStringArg("-chr");
		int start = p.getIntArg("-s");
		int end = p.getIntArg("-e");
		boolean useFragments = p.getBooleanArg("-pe");
		
		TranscriptomeSpace t = new TranscriptomeSpace(BEDFileParser.loadDataByChr(new File(bedFile)));
		GenomicSpace g = new GenomicSpace(chrSizeFile);
		
		AlignmentModel transcriptomeAlignments = new AlignmentModel(bamFile, t, useFragments);
		AlignmentModel genomeAlignments = new AlignmentModel(bamFile, g, useFragments);
		
		Gene testGene = new Gene(chrName, start, end);
		
		double tCount = transcriptomeAlignments.getCount(testGene,true);
		double gCount = genomeAlignments.getCount(testGene,true);
		
		System.out.println("Window: " + chrName + ":" + start + "-" + end);
		System.out.println("Count in transcriptome space: " + tCount);
		System.out.println("Count in genome space: " + gCount);
		
	}

}
