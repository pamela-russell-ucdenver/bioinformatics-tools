/**
 * 
 */
package util.programs.bed;

import guttmanlab.core.util.CommandLineParser;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.log4j.Logger;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Gene;

import broad.pda.annotation.BEDFileParser;

/**
 * @author prussell
 *
 */
public class BedFileFilter {

	Map<String, Collection<Gene>> regions;
	private static Logger logger = Logger.getLogger(BedFileFilter.class.getName());
	
	/**
	 * @param mainGenesBed Main genes to filter
	 * @throws IOException
	 */
	public BedFileFilter(String mainGenesBed) throws IOException {
		logger.info("");
		logger.info("Loading genes from file " + mainGenesBed + ".");
		regions = BEDFileParser.loadDataByChr(new File(mainGenesBed));
		int numRegions = 0;
		for(String chr : regions.keySet()) {
			numRegions += regions.get(chr).size();
		}
		logger.info("Loaded " + numRegions + " genes.");
	}
	
	/**
	 * Write genes with names on the list to file
	 * @param geneNamesToKeepFile File containing list of gene names
	 * @param outBedFile Output bed file
	 * @throws IOException
	 */
	public void filterGenesByName(String geneNamesToKeepFile, String outBedFile) throws IOException {

		logger.info("");
		logger.info("Filtering genes by name. Keeping genes in file " + geneNamesToKeepFile + ".");
		
        FileReader genereader = new FileReader(new File(geneNamesToKeepFile));
        BufferedReader genebuffered = new BufferedReader(genereader);
        TreeSet<String> geneNames = new TreeSet<String>();

        while(genebuffered.ready()) {
             String line = genebuffered.readLine();
             geneNames.add(line);
        }

        FileWriter writer = new FileWriter(outBedFile);

        for(String chr : regions.keySet()) {
        	for(Gene gene : regions.get(chr)) {
        		if(geneNames.contains(gene.getName())) {
        			writer.write(gene.toBED() + "\n");
        		}
        	}
        }
        
        writer.close();
        genereader.close();
        genebuffered.close();
        
        logger.info("Done filtering genes by name. Wrote to " + outBedFile + ".");

	}
	
	/**
	 * Filter genes by overlap with another set of genes and write filtered genes to file
	 * @param otherGenesBedFile The other set of genes
	 * @param outBedFile Output bed file
	 * @param keepOverlappers Keep overlappers and remove non-overlappers (if false, keep non-overlappers and remove overlappers)
	 * @param ignoreOrientation Ignore orientation for overlap. If false, genes must have same or unknown orientation to be called overlap.
	 * @param ignoreIfSelf Do not consider overlap if gene is equal to other gene
	 * @throws IOException
	 */
	public void filterGenesByOverlap(String otherGenesBedFile, String outBedFile, boolean keepOverlappers, boolean ignoreOrientation, boolean ignoreIfSelf) throws IOException {
		
		logger.info("");
		logger.info("Filtering genes by overlap with genes in file " + otherGenesBedFile + ".");
		if(keepOverlappers) logger.info("Keeping overlappers.");
		else logger.info("Removing overlappers.");
		
		FileWriter w = new FileWriter(outBedFile);
		Map<String, Collection<Gene>> otherGenes = BEDFileParser.loadDataByChr(new File(otherGenesBedFile));
		Map<String, Collection<Gene>> filteredGenes = filterGenesByOverlap(regions, otherGenes, keepOverlappers, ignoreOrientation, ignoreIfSelf);
		for(String chr : filteredGenes.keySet()) {
			for(Gene gene : filteredGenes.get(chr)) {
				w.write(gene.toBED() + "\n");
			}
		}
		w.close();
		
		logger.info("Done filtering genes by overlap. Wrote results to " + outBedFile + ".");
		
	}
	
	/**
	 * Filter genes by overlap with another set of genes
	 * @param genes The set of genes to filter (by chromosome name)
	 * @param otherGenes The other set of genes (by chromosome name)
	 * @param keepOverlappers Keep overlappers and remove non-overlappers (if false, keep non-overlappers and remove overlappers)
	 * @param ignoreOrientation Ignore orientation for overlap. If false, genes must have same or unknown orientation to be called overlap.
	 * @param ignoreIfSelf Do not consider overlap if gene is equal to other gene
	 * @return The filtered set of genes
	 */
	public static Map<String, Collection<Gene>> filterGenesByOverlap(Map<String, Collection<Gene>> genes, Map<String, Collection<Gene>> otherGenes, boolean keepOverlappers, boolean ignoreOrientation, boolean ignoreIfSelf) {
		Map<String, Collection<Gene>> rtrn = new TreeMap<String, Collection<Gene>>();
		for(String chr : genes.keySet()) {
			logger.info(chr);
			Collection<Gene> chrOtherGenes = otherGenes.get(chr);
			Collection<Gene> chrGenes = new TreeSet<Gene>();
			chrGenes.addAll(genes.get(chr));
			for(Gene gene : genes.get(chr)) {
				boolean foundOverlap = false;
				for(Gene otherGene : chrOtherGenes) {
					if(gene.getEnd() < otherGene.getStart() || gene.getStart() > otherGene.getEnd()) {
						continue;
					}
					boolean ignore = ignoreIfSelf && gene.equals(otherGene);
					boolean overlaps = gene.overlaps(otherGene, ignoreOrientation) && !ignore;
					if(overlaps) foundOverlap = true;
				}
				boolean keep = (foundOverlap && keepOverlappers) || (!foundOverlap && !keepOverlappers);
				if(!keep) {
					chrGenes.remove(gene);
				}
			}
			rtrn.put(chr, chrGenes);
		}		
		return rtrn;
	}
	
	/**
	 * Subtract segments from the genes and write result to bed file
	 * @param otherGenesBedFile Bed file of segments to subtract
	 * @param outBedFile Output bed file
	 * @param separateBlocks After subtracting other genes from an annotation, leave the result as separate single-exon annotations. If false, keep the block structure in a single annotation.
	 * @throws IOException
	 */
	public void subtractSegmentsFromGenes(String otherGenesBedFile, String outBedFile, boolean separateBlocks) throws IOException {
		
		logger.info("");
		logger.info("Subtracting segments in file " + otherGenesBedFile + ".");
		
		Map<String, Collection<Gene>> otherGenes = BEDFileParser.loadDataByChr(new File(otherGenesBedFile));
		Collection<Annotation> subtracted = subtractSegmentsFromGenes(regions, otherGenes);
		FileWriter w = new FileWriter(outBedFile);
		for(Annotation region : subtracted) {
			if(separateBlocks) {
				Collection<? extends Annotation> exons = region.getBlocks();
				int blockNum = 0;
				for(Annotation exon : exons) {
					exon.setName(region.getName() + "_block_" + blockNum);
					w.write(exon.toBED() + "\n");
					blockNum++;
				}
			} else {
				w.write(region.toBED() + "\n");
			}
		}
		w.close();
		
		logger.info("Done subtracting segments. Wrote results to file " + outBedFile + ".");
		
	}
	
	/**
	 * Subtract overlapping segments from genes
	 * Keep the rest of the gene model
	 * @param genes The genes
	 * @param otherGenes The segments to subtract
	 * @return The genes minus the segments
	 */
	public static Collection<Annotation> subtractSegmentsFromGenes(Map<String, Collection<Gene>> genes, Map<String, Collection<Gene>> otherGenes) {
		Collection<Annotation> rtrn = new TreeSet<Annotation>();
		for(String chr : genes.keySet()) {
			logger.info(chr);
			if(!otherGenes.containsKey(chr)) {
				rtrn.addAll(genes.get(chr));
				continue;
			}
			if(otherGenes.get(chr).isEmpty()) {
				rtrn.addAll(genes.get(chr));
				continue;
			}
			Collection<Gene> chrOtherGenes = otherGenes.get(chr);
			for(Gene gene : genes.get(chr)) {
				Annotation newGene = gene.minus(chrOtherGenes);
				if(newGene != null) {
					if(newGene.size() < 1) continue;
					newGene.setName(gene.getName());
					rtrn.add(newGene);
				}
			}
		}
		return rtrn;
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		
        CommandLineParser p = new CommandLineParser();
        p.addStringArg("-b", "Bed file of genes to filter", true);
        p.addStringArg("-gk", "File containing list of gene names to keep", false, null);
        p.addStringArg("-og", "Output bed file for genes filtered by gene name. Requires -gk.",false, null);
        p.addStringArg("-go", "Bed file of other genes for overlap filter", false, null);
        p.addStringArg("-oo", "Output bed file for whole genes filtered by overlap. Requires -go.", false, null);
        p.addStringArg("-os", "Output bed file for genes with other segments subtracted from gene model. Requires -go.", false, null);
        p.addBooleanArg("-ko", "For overlap filter, keep overlappers and filter non-overlappers", false, false);
        p.addBooleanArg("-io", "Ignore orientation for overlap filter", false, false);
        p.addBooleanArg("-is", "Ignore gene for overlap feature if gene is equal to other gene", false, false);
        p.addBooleanArg("-ss", "For subtraction of overlappers, separate resulting annotations into contiguous blocks", false, false);
        p.parse(args);
        String mainGenesBed = p.getStringArg("-b");
        String geneNamesToKeepFile = p.getStringArg("-gk");
        String outBedFile_filterByGeneName = p.getStringArg("-og");
        String overlapBedFile = p.getStringArg("-go");
        String outBedFile_filterByOverlap = p.getStringArg("-oo");
        String outBedFile_subtract = p.getStringArg("-os");
        boolean keepOverlappers = p.getBooleanArg("-ko");
        boolean ignoreOrientation = p.getBooleanArg("-io");
        boolean ignoreIfSelf = p.getBooleanArg("-is");
        boolean separate = p.getBooleanArg("-ss");
        
        // Check for valid command line
        if(outBedFile_filterByGeneName != null) {
        	if(geneNamesToKeepFile == null) {
        		throw new IllegalArgumentException("Must provide -gk.");
        	}
        }
        if(outBedFile_filterByOverlap != null || outBedFile_subtract != null) {
        	if(overlapBedFile == null) {
        		throw new IllegalArgumentException("Must provide -go.");
        	}
        }
        
        // Instantiate
        BedFileFilter bff = new BedFileFilter(mainGenesBed);
        
        // Filter genes by name
        if(outBedFile_filterByGeneName != null) {
        	bff.filterGenesByName(geneNamesToKeepFile, outBedFile_filterByGeneName);
        }
        	
        // Filter genes by overlap
        if(outBedFile_filterByOverlap != null) {
         	bff.filterGenesByOverlap(overlapBedFile, outBedFile_filterByOverlap, keepOverlappers, ignoreOrientation, ignoreIfSelf);
        }
        
        // Subtract segments from genes
        if(outBedFile_subtract != null) {
        	bff.subtractSegmentsFromGenes(overlapBedFile, outBedFile_subtract, separate);
        }
        
 	}

}
