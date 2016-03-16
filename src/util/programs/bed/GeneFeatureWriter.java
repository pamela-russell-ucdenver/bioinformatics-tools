/**
 * 
 */
package util.programs.bed;

import guttmanlab.core.util.CommandLineParser;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.TreeSet;

import org.apache.log4j.Logger;

import broad.pda.annotation.BEDFileParser;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.BasicAnnotation;
import nextgen.core.annotation.Gene;
import nextgen.core.annotation.Annotation.Strand;
import nextgen.core.utils.AnnotationUtils;

/**
 * @author prussell
 *
 */
public class GeneFeatureWriter {

	private Collection<Gene> genes;
	private static Logger logger = Logger.getLogger(GeneFeatureWriter.class.getName());
	
	private GeneFeatureWriter(String bedFile) throws IOException {
		genes = BEDFileParser.loadData(new File(bedFile));
	}
	
	private Collection<Gene> getOriginalGenes() {
		return genes;
	}
	
	private Collection<Gene> getIndividualIntrons() {
		logger.info("Getting all individual introns...");
		Collection<Gene> rtrn = new TreeSet<Gene>();
		for(Gene gene : genes) {
			Collection<? extends Annotation> introns = gene.getIntronSet();
			for(Annotation intron : introns) {
				Gene intronGene = new Gene(intron);
				rtrn.add(intronGene);
			}
		}
		logger.info("Got " + rtrn.size() + " introns.");
		return rtrn;
	}
	
	private Collection<Gene> getIndividualExons() {
		logger.info("Getting all individual exons...");
		Collection<Gene> rtrn = new TreeSet<Gene>();
		for(Gene gene : genes) {
			int exonNumber = 0;
			Collection<? extends Annotation> exons = gene.getExonSet();
			for(Annotation exon : exons) {
				Gene exonGene = new Gene(exon);
				exonGene.setName(gene.getName() + "_exon_" + exonNumber);
				exonNumber++;
				rtrn.add(exonGene);
			}
		}
		logger.info("Got " + rtrn.size() + " exons.");
		return rtrn;
	}
	
	private Collection<Gene> getNonCodingRNAs() {
		logger.info("Getting all non-coding genes...");
		Collection<Gene> rtrn = new TreeSet<Gene>();
		for(Gene gene : genes) {
			if(!gene.isCodingGene()) {
				Gene nc = gene.copy();
				nc.setName(gene.getName() + "_noncoding");
				rtrn.add(nc);
			}
		}
		logger.info("Got " + rtrn.size() + " ncRNAs.");
		return rtrn;
	}
	
	private Collection<Gene> getAllIntronSets() {
		logger.info("Getting sets of introns for all genes...");
		Collection<Gene> rtrn = new TreeSet<Gene>();
		for(Gene gene : genes) {
			Gene introns = gene.getIntrons();
			if(introns != null) {
				introns.setName(gene.getName() + "_all_introns");
				introns.setCDSStart(introns.getStart());
				introns.setCDSEnd(introns.getStart());
				rtrn.add(introns);
			}
		}
		logger.info("Got " + rtrn.size() + " intron sets.");
		return rtrn;
	}
	
	private Collection<Gene> getAllCDSs() {
		logger.info("Getting all CDSs...");
		Collection<Gene> rtrn = new TreeSet<Gene>();
		for(Gene gene : genes) {
			if(!gene.isCodingGene()) {
				continue;
			}
			Gene cds = gene.getCDS();
			if(cds == null) {
				continue;
			}
			cds.setName(gene.getName() + "_CDS");
			rtrn.add(cds);
		}
		logger.info("Got " + rtrn.size() + " CDSs.");
		return rtrn;
	}
	
	
	private Collection<Gene> getAll3UTRs() {
		logger.info("Getting all 3'-UTRs...");
		Collection<Gene> rtrn = new TreeSet<Gene>();
		for(Gene gene : genes) {
			if(!gene.isCodingGene()) {
				continue;
			}
			Gene utr3 = gene.get3UTRGene();
			if(utr3 == null) {
				continue;
			}
			utr3.setName(gene.getName() + "_3UTR");
			rtrn.add(utr3);
		}
		logger.info("Got " + rtrn.size() + " 3'-UTRs.");
		return rtrn;
	}
	
	private Collection<Gene> getAll5UTRs() {
		logger.info("Getting all 5'-UTRs...");
		Collection<Gene> rtrn = new TreeSet<Gene>();
		for(Gene gene : genes) {
			if(!gene.isCodingGene()) {
				continue;
			}
			Gene utr5 = gene.get5UTRGene();
			if(utr5 == null) {
				continue;
			}
			utr5.setName(gene.getName() + "_5UTR");
			rtrn.add(utr5);
		}
		logger.info("Got " + rtrn.size() + " 5'-UTRs.");
		return rtrn;
	}
	
	/**
	 * Get the "beginning" of each gene
	 * If window size is larger than first exon, return entire first exon only
	 * @param windowSize Window size to get
	 * @return Window beginning at TSS and extending into the first exon
	 */
	private Collection<Gene> getWindowDownstreamOfEachTSS(int windowSize) {
		Collection<Gene> rtrn = new TreeSet<Gene>();
		for(Gene gene : genes) {
			Strand orientation = gene.getOrientation();
			if(orientation.equals(Strand.UNKNOWN)) {
				logger.warn("Skipping gene " + gene.getName() + " because orientation is unknown.");
			}
			Iterator<? extends Annotation> iter = null;
			TreeSet<Annotation> exons = new TreeSet<Annotation>();
			exons.addAll(gene.getExonSet());
			boolean plusStrand = orientation.equals(Strand.POSITIVE);
			if(plusStrand) {
				iter = exons.iterator();
			} else {
				iter = exons.descendingIterator();
			}
			Annotation exon = iter.next();
			if(exon.getSize() <= windowSize) {
				Gene window = new Gene(exon);
				window.setName(gene.getName() + "_first_" + windowSize + "bp_downstream_of_TSS");
				rtrn.add(window);
				continue;
			}
			int windowLastPos = plusStrand ? exon.getStart() + windowSize - 1 : exon.getEnd() - windowSize;
			int windowFirstPos = plusStrand ? exon.getStart() : exon.getEnd() - 1;
			int windowStart = Math.min(windowFirstPos, windowLastPos);
			int windowEnd = Math.max(windowFirstPos, windowLastPos) + 1;
			Gene window = new Gene(new BasicAnnotation(gene.getChr(), windowStart, windowEnd, orientation));
			window.setName(gene.getName() + "_first_" + windowSize + "bp_downstream_of_TSS");
			rtrn.add(window);
		}
		return rtrn;
	}
	
	/**
	 * Get a window upstream of each exon exon junction in all genes
	 * @param windowSize Window size
	 * @param distanceUpstreamOfExon3PrimeEnd Distance upstream of 3' end of exon (number of bases to leave off of 3' end)
	 * @return Windows for all exon exon junctions of all genes
	 */
	private Collection<Gene> getWindowUpstreamOfEachExonExonJunction(int windowSize, int distanceUpstreamOfExon3PrimeEnd) {
		Collection<Gene> rtrn = new TreeSet<Gene>();
		for(Gene gene : genes) {
			if(gene.numBlocks() < 2) continue;
			Strand orientation = gene.getOrientation();
			if(orientation.equals(Strand.UNKNOWN)) {
				logger.warn("Skipping gene " + gene.getName() + " because orientation is unknown.");
			}
			Iterator<? extends Annotation> iter = null;
			TreeSet<Annotation> exons = new TreeSet<Annotation>();
			exons.addAll(gene.getExonSet());
			boolean plusStrand = orientation.equals(Strand.POSITIVE);
			if(plusStrand) {
				iter = exons.iterator();
			} else {
				iter = exons.descendingIterator();
			}
			int exonNumber = 1;
			while(iter.hasNext()) {
				Annotation exon = iter.next();
				if(!iter.hasNext()) {
					// Don't use the last exon
					continue;
				}
				if(exon.getSize() < windowSize + distanceUpstreamOfExon3PrimeEnd) {
					continue;
				}
				int lastExonPos = plusStrand ? exon.getEnd() - 1 : exon.getStart();
				int windowLastPos = plusStrand ? lastExonPos - distanceUpstreamOfExon3PrimeEnd : lastExonPos + distanceUpstreamOfExon3PrimeEnd;
				int windowFirstPos = plusStrand ? windowLastPos - windowSize + 1 : windowLastPos + windowSize - 1;
				int windowStart = Math.min(windowFirstPos, windowLastPos);
				int windowEnd = Math.max(windowFirstPos, windowLastPos) + 1;
				Gene window = new Gene(new BasicAnnotation(gene.getChr(), windowStart, windowEnd, orientation));
				window.setName(gene.getName() + "_exon_exon_junction_" + exonNumber);
				rtrn.add(window);
				exonNumber++;
			}
		}
		return rtrn;
	}
	
	/**
	 * Write individual exons to bed file
	 * @param inputBed Input bed
	 * @param outputBed Output bed file of exons
	 * @throws IOException
	 */
	public static void writeIndividualExons(String inputBed, String outputBed) throws IOException {
		GeneFeatureWriter w = new GeneFeatureWriter(inputBed);
		Collection<Gene> exons = w.getIndividualExons();
		Collection<Annotation> e = new TreeSet<Annotation>();
		e.addAll(exons);
		w.writeAll(e, outputBed);
	}
	
	private void writeAll(Collection<Annotation> features, String fileName) throws IOException {
		logger.info("Writing " + features.size() + " genes to file " + fileName);
		FileWriter w = new FileWriter(fileName);
		for(Annotation feature : features) {
			w.write(feature.toBED() + "\n");
		}
		w.close();
		logger.info("Done writing file.");
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-i", "Input bed file", true);
		p.addBooleanArg("--orig", "Write original genes", false, false);
		p.addBooleanArg("--intron_sets", "Write intron set for each gene", false, false);
		p.addBooleanArg("--ind_introns", "Write each individual intron for each gene", false, false);
		p.addBooleanArg("--ind_exons", "Write each individual exon for each gene", false, false);
		p.addBooleanArg("--utr3", "Write 3'-UTRs", false, false);
		p.addBooleanArg("--utr5", "Write 5'-UTRs", false, false);
		p.addBooleanArg("--cds", "Write CDSs", false, false);
		p.addBooleanArg("--nc", "Write original non-coding genes", false, false);
		p.addBooleanArg("--eej", "Write window upstream of each exon exon junction", false, false);
		p.addIntArg("--eejw", "Window size for -eej", false, 14);
		p.addIntArg("--eejd", "Distance upstream of exon exon junction for --eej", false, 17);
		p.addBooleanArg("--dtss", "Write window downstream of each TSS", false, false);
		p.addIntArg("--dtssw", "Window size for --dtss", false, 200);
		p.addBooleanArg("-mo", "Merge overlappers", false, false);
		p.addStringArg("-o", "Output bed file", true);
		p.parse(args);
		
		String inputBed = p.getStringArg("-i");
		GeneFeatureWriter gfw = new GeneFeatureWriter(inputBed);
		String outputBed = p.getStringArg("-o");
		boolean intronSets = p.getBooleanArg("--intron_sets");
		boolean indIntrons = p.getBooleanArg("--ind_introns");
		boolean indExons = p.getBooleanArg("--ind_exons");
		boolean orig = p.getBooleanArg("--orig");
		boolean utr3 = p.getBooleanArg("--utr3");
		boolean cds = p.getBooleanArg("--cds");
		boolean utr5 = p.getBooleanArg("--utr5");
		boolean nc = p.getBooleanArg("--nc");
		boolean eej = p.getBooleanArg("--eej");
		int eejw = p.getIntArg("--eejw");
		int eejd = p.getIntArg("--eejd");
		boolean dtss = p.getBooleanArg("--dtss");
		int dtssw = p.getIntArg("--dtssw");
		boolean merge = p.getBooleanArg("-mo");
		
		Collection<Annotation> genesToWrite = new TreeSet<Annotation>();
		if(orig) genesToWrite.addAll(gfw.getOriginalGenes());
		if(intronSets) genesToWrite.addAll(gfw.getAllIntronSets());
		if(indIntrons) genesToWrite.addAll(gfw.getIndividualIntrons());
		if(indExons) genesToWrite.addAll(gfw.getIndividualExons());
		if(utr3) genesToWrite.addAll(gfw.getAll3UTRs());
		if(utr5) genesToWrite.addAll(gfw.getAll5UTRs());
		if(nc) genesToWrite.addAll(gfw.getNonCodingRNAs());
		if(eej) genesToWrite.addAll(gfw.getWindowUpstreamOfEachExonExonJunction(eejw, eejd));
		if(dtss) genesToWrite.addAll(gfw.getWindowDownstreamOfEachTSS(dtssw));
		if(cds) genesToWrite.addAll(gfw.getAllCDSs());
		if(merge) {
			TreeSet<Annotation> t = new TreeSet<Annotation>();
			t.addAll(genesToWrite);
			Collection<Annotation> c = AnnotationUtils.mergeOverlappingBlocks(t);
			gfw.writeAll(c, outputBed);
		} else {
			gfw.writeAll(genesToWrite, outputBed);
		}
		
	}

}
