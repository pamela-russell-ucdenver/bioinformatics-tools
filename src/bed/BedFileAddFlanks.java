package bed;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.BlockedAnnotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.AnnotationCollection;
import guttmanlab.core.util.CommandLineParser;

public class BedFileAddFlanks {
	
	private static Annotation addFlanks(Annotation region, int flankSize) {
		SingleInterval firstBlock = null;
		SingleInterval lastBlock = null;
		int firstBlockStart = Integer.MAX_VALUE;
		int lastBlockEnd = Integer.MIN_VALUE;
		Iterator<SingleInterval> iter = region.getBlocks();
		while(iter.hasNext()) {
			SingleInterval block = iter.next();
			int start = block.getReferenceStartPosition();
			int end = block.getReferenceEndPosition();
			if(start < firstBlockStart) {
				firstBlock = block;
				firstBlockStart = start;
			}
			if(end > lastBlockEnd) {
				lastBlock = block;
				lastBlockEnd = end;
			}
		}
		SingleInterval firstBlockExpanded = new SingleInterval(firstBlock.getReferenceName(), firstBlock.getReferenceStartPosition() - flankSize, firstBlock.getReferenceEndPosition(), firstBlock.getOrientation());
		SingleInterval lastBlockExpanded = new SingleInterval(lastBlock.getReferenceName(), lastBlock.getReferenceStartPosition(), lastBlock.getReferenceEndPosition() + flankSize, lastBlock.getOrientation());
		Collection<Annotation> newBlocks = new ArrayList<Annotation>();
		Iterator<SingleInterval> iter2 = region.getBlocks();
		while(iter2.hasNext()) {
			SingleInterval block = iter2.next();
			if(block.equals(firstBlock) && block.equals(lastBlock)) {
				SingleInterval newBlock = new SingleInterval(firstBlock.getReferenceName(), firstBlockExpanded.getReferenceStartPosition(), lastBlockExpanded.getReferenceEndPosition(), firstBlock.getOrientation());
				newBlocks.add(newBlock);
				continue;
			}
			if(block.equals(firstBlock)) {
				newBlocks.add(firstBlockExpanded);
				continue;
			}
			if(block.equals(lastBlock)) {
				newBlocks.add(lastBlockExpanded);
				continue;
			}
			newBlocks.add(block);
		}
		String newName = region.getName() + "_plus_" + flankSize + "bp_flanks";
		return new BlockedAnnotation(newBlocks, newName);
	}
	
	private static void addFlanksAndWrite(String inputBed, String chrSizeFile, int flankSize, String outputBed) throws IOException {
		AnnotationCollection<Gene> genes = BEDFileIO.loadFromFile(inputBed, chrSizeFile);
		FileWriter w = new FileWriter(outputBed);
		Iterator<Gene> iter = genes.sortedIterator();
		while(iter.hasNext()) {
			Gene gene = iter.next();
			Annotation expanded = addFlanks(gene, flankSize);
			w.write(expanded.toBED() + "\n");
		}
		w.close();
	}
	
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-i", "Input bed", true);
		p.addStringArg("-o", "Output bed", true);
		p.addStringArg("-c", "Chr size file", true);
		p.addIntArg("-f", "Flank size to add", true);
		p.parse(args);
		String inputBed = p.getStringArg("-i");
		String outputBed = p.getStringArg("-o");
		String chrSizeFile = p.getStringArg("-c");
		int flankSize = p.getIntArg("-f");
		
		addFlanksAndWrite(inputBed, chrSizeFile, flankSize, outputBed);

	}

}
