package annotation;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.core.parser.CommandLineParser;
import broad.pda.annotation.BEDFileParser;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Annotation.Strand;
import nextgen.core.annotation.Gene;


public class ConstitutiveExonFinder {
	
	private Map<String, Collection<Gene>> genes;
	private Map<String, Collection<Annotation>> isoformsByGeneName; // Isoforms are grouped if they have the same name
	private Map<String, Collection<Annotation>> exonsByGeneName;
	private Map<String, Collection<Annotation>> constitutiveExonsByGeneName;
	
	/**
	 * @param geneBed Genes to use for constitutive exons
	 * @throws IOException
	 */
	private ConstitutiveExonFinder(String geneBed) throws IOException {
		genes = BEDFileParser.loadDataByChr(geneBed);
		populateIsoformsByGeneName();
		populateExonsByGeneName();
		findConstitutiveExons();
	}
	
	/**
	 * Write constitutive exons to a file
	 * @param outBed Output bed file
	 * @throws IOException
	 */
	private void writeConstitutiveExons(String outBed) throws IOException {
		FileWriter w = new FileWriter(outBed);
		for(String gene : constitutiveExonsByGeneName.keySet()) {
			for(Annotation exon : constitutiveExonsByGeneName.get(gene)) {
				w.write(exon.toBED() + "\n");
			}
		}
		w.close();
	}
	
	/**
	 * Group transcripts that have the same name
	 */
	private void populateIsoformsByGeneName() {
		isoformsByGeneName = new TreeMap<String, Collection<Annotation>>();
		for(String chr : genes.keySet()) {
			for(Gene gene : genes.get(chr)) {
				String geneName = gene.getName();
				if(!isoformsByGeneName.containsKey(geneName)) {
					isoformsByGeneName.put(geneName, new TreeSet<Annotation>());
				}
				isoformsByGeneName.get(geneName).add(gene);
			}
		}
	}
	
	/**
	 * Reset the name and erase info that might affect equals() but doesn't matter here
	 * @param block Annotation to modify
	 * @param newName New name for annotation
	 */
	private void eraseExtraInfo(Annotation block, String newName) {
		block.setName(newName);
		block.setScore(0);
	}
	
	/**
	 * Transcripts are grouped if they have the same name
	 */
	private void populateExonsByGeneName() {
		exonsByGeneName = new TreeMap<String, Collection<Annotation>>();
		for(String chr : genes.keySet()) {
			for(Gene gene : genes.get(chr)) {
				String geneName = gene.getName();
				if(!exonsByGeneName.containsKey(geneName)) {
					exonsByGeneName.put(geneName, new TreeSet<Annotation>());
				}
				List<? extends Annotation> blocks = gene.getBlocks();
				for(Annotation block : blocks) {
					eraseExtraInfo(block, gene.getName());
					exonsByGeneName.get(geneName).add(block);
				}
			}
		}
	}
	
	/**
	 * Populate the set of constitutive exons
	 */
	private void findConstitutiveExons() {
		constitutiveExonsByGeneName = new TreeMap<String, Collection<Annotation>>();
		for(String geneName : exonsByGeneName.keySet()) {
			if(!constitutiveExonsByGeneName.containsKey(geneName)) {
				constitutiveExonsByGeneName.put(geneName, new TreeSet<Annotation>());
			}
			for(Annotation exon : exonsByGeneName.get(geneName)) {
				if(exonIsConstitutive(geneName, exon)) {
					constitutiveExonsByGeneName.get(geneName).add(exon);
				}
			}
		}
	}
	
	/**
	 * Check if the two annotations are equivalent single blocks
	 * @param a1
	 * @param a2
	 * @return True iff the two annotations each have a single block with same chr, start, end, strand
	 */
	private boolean singleIntervalEquals(Annotation a1, Annotation a2) {
		if(a1.numBlocks() != 1 || a2.numBlocks() != 1) {
			throw new IllegalArgumentException("Both annotations must have one block");
		}
		String chr1 = a1.getChr();
		String chr2 = a2.getChr();
		if(!chr1.equals(chr2)) return false;
		int start1 = a1.getStart();
		int start2 = a2.getStart();
		if(start1 != start2) return false;
		int end1 = a1.getEnd();
		int end2 = a2.getEnd();
		if(end1 != end2) return false;
		Strand strand1 = a1.getOrientation();
		Strand strand2 = a2.getOrientation();
		if(!strand1.equals(strand2)) return false;
		return true;
	}
	
	/**
	 * Check if an exon is shared by all isoforms of the specified gene
	 * @param geneName 
	 * @param exon
	 * @return
	 */
	private boolean exonIsConstitutive(String geneName, Annotation exon) {
		for(Annotation isoform : isoformsByGeneName.get(geneName)) {
			boolean hasExon = false;
			Iterator<? extends Annotation> blockIter = isoform.getBlocks().iterator();
			while(blockIter.hasNext()) {
				Annotation block = blockIter.next();
				if(singleIntervalEquals(exon, block)) {
					hasExon = true;
					break;
				}
			}
			if(!hasExon) return false;
		}
		return true;
	}
	
	
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-b", "Bed file of genes", true);
		p.addStringArg("-o", "Output bed file of constitutive exons", true);
		p.parse(args);
		String bed = p.getStringArg("-b");
		String out = p.getStringArg("-o");
		
		ConstitutiveExonFinder finder = new ConstitutiveExonFinder(bed);
		finder.writeConstitutiveExons(out);

	}

}
