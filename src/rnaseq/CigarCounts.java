/**
 * 
 */
package rnaseq;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;

import org.apache.log4j.Logger;

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.TextCigarCodec;
import net.sf.samtools.util.CloseableIterator;
import nextgen.core.alignment.Alignment;
import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Gene;
import nextgen.core.coordinatesystem.GenomicSpace;
import nextgen.core.coordinatesystem.TranscriptomeSpace;
import nextgen.core.model.AlignmentModel;


import broad.core.parser.CommandLineParser;
import broad.pda.annotation.BEDFileParser;

/**
 * @author prussell
 *
 */
public class CigarCounts {

	private AlignmentModel transcriptomeData;
	private AlignmentModel genomicData;
	private Map<String, Collection<Gene>> genesByChr;
	private Collection<Gene> genes;
	private Map<String, Collection<Gene>> featuresByChr;
	private Collection<Gene> features;
	private static Logger logger = Logger.getLogger(CigarCounts.class.getName());
	
	private CigarCounts(String bamFile, String chrSizeFile, String bedGenes, String bedFeatures) throws IOException {
		genesByChr = BEDFileParser.loadDataByChr(new File(bedGenes));
		genes = new TreeSet<Gene>();
		for(String chr : genesByChr.keySet()) {
			genes.addAll(genesByChr.get(chr));
		}
		featuresByChr = BEDFileParser.loadDataByChr(new File(bedFeatures));
		features = new TreeSet<Gene>();
		for(String chr : featuresByChr.keySet()) {
			features.addAll(featuresByChr.get(chr));
		}
		transcriptomeData = new AlignmentModel(bamFile, new TranscriptomeSpace(genesByChr), false);
		genomicData = new AlignmentModel(bamFile, new GenomicSpace(chrSizeFile), false);
	}
	
	private static boolean cigarHasDeletion(Alignment align) {
		String cigarString = align.toSAMRecord().getCigarString();
		Cigar cigar = TextCigarCodec.getSingleton().decode(cigarString);
		List<CigarElement> cigarElements = cigar.getCigarElements();
		for(CigarElement cigarElement : cigarElements) {
			if(cigarElement.getOperator().equals(CigarOperator.DELETION)) {
				return true;
			}
		}
		return false;
	}
	
	private double pctOverlappersWithDeletion() {
		logger.info("Calculating percentage of reads overlapping genomic space that contain a deletion.");
		double hasDeletion = 0;
		double noDeletion = 0;
		int numDone = 0;
		CloseableIterator<Alignment> iter = genomicData.getReadIterator();
		while(iter.hasNext()) {
			numDone++;
			if(numDone % 10000000 == 0) {
				logger.info("Finished " + numDone + " reads. Percentage with deletion: " + hasDeletion / (hasDeletion + noDeletion));
			}
			Alignment align = iter.next();
			if(cigarHasDeletion(align)) {
				hasDeletion += align.getWeight();
			} else {
				noDeletion += align.getWeight();
			}
		}
		logger.info("Has deletion: " + hasDeletion);
		logger.info("No deletion: " + noDeletion);
		double pct = hasDeletion / (hasDeletion + noDeletion);
		logger.info("Percentage with deletion: " + pct );
		return pct;
	}
	
	@SuppressWarnings("unused")
	private double pctOverlappersWithDeletion(Collection<Gene> regions) {
		return pctOverlappersWithDeletion(regions, false);
	}
	
	private double pctOverlappersWithDeletion(Collection<Gene> regions, boolean fullyContained) {
		logger.info("Calculating percentage of reads overlapping " + regions.size() + " regions that contain a deletion.");
		double hasDeletion = 0;
		double noDeletion = 0;
		int numDone = 0;
		for(Annotation region : regions) {
			numDone++;
			if(numDone % 1000 == 0) {
				logger.info("Finished " + numDone + " regions. Percentage with deletion: " + hasDeletion / (hasDeletion + noDeletion));
			}
			CloseableIterator<Alignment> iter = transcriptomeData.getOverlappingReads(region, fullyContained);
			while(iter.hasNext()) {
				Alignment align = iter.next();
				if(cigarHasDeletion(align)) {
					hasDeletion += align.getWeight();
				} else {
					noDeletion += align.getWeight();
				}
			}
		}
		logger.info("Has deletion: " + hasDeletion);
		logger.info("No deletion: " + noDeletion);
		double pct = hasDeletion / (hasDeletion + noDeletion);
		logger.info("Percentage with deletion: " + pct);
		return pct;
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-b", "Bam file", true);
		p.addStringArg("-g", "Gene bed file", true);
		p.addStringArg("-f", "Feature bed file", true);
		p.addStringArg("-c", "Chromosome size file", true);
		p.parse(args);
		String bamFile = p.getStringArg("-b");
		String geneBedFile = p.getStringArg("-g");
		String featureBedFile = p.getStringArg("-f");
		String chrSizeFile = p.getStringArg("-c");

		CigarCounts c = new CigarCounts(bamFile, chrSizeFile, geneBedFile, featureBedFile);
		logger.info("Calculating percentage of reads overlaping features that contain a deletion...");
		@SuppressWarnings("unused")
		double pctWithDeletionFeatures = c.pctOverlappersWithDeletion(c.features, false);
		logger.info("Calculating percentage of reads overlapping genes that contain a deletion...");
		@SuppressWarnings("unused")
		double pctWithDeletionGenes = c.pctOverlappersWithDeletion(c.genes, false);
		logger.info("Calculating percentage of all reads containing a deletion...");
		@SuppressWarnings("unused")
		double pctWithDeletionAll = c.pctOverlappersWithDeletion();
		
	}

}
