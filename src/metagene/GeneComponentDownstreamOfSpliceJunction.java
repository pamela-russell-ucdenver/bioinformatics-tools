package metagene;

import java.util.Collection;
import java.util.Iterator;
import java.util.TreeSet;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.BasicAnnotation;
import nextgen.core.annotation.Gene;
import nextgen.core.annotation.Annotation.Strand;

/**
 * @author prussell
 *
 */
public class GeneComponentDownstreamOfSpliceJunction extends AbstractGeneComponent {

	private int numBases;
	
	/**
	 * @param numBasesUpstream Number of bases to include
	 */
	public GeneComponentDownstreamOfSpliceJunction(int numBasesUpstream) {
		numBases = numBasesUpstream;
	}
	
	@Override
	public String getName() {
		return "downstream_of_splice_junction_" + numBases;
	}

	@Override
	public Collection<Gene> getComponent(Gene gene) {
	
		if(gene.numBlocks() < 2) {
			return null;
		}
		Strand orientation = gene.getOrientation();
		if(orientation.equals(Strand.UNKNOWN)) {
			return null;
		}
		
		Collection<Gene> rtrn = new TreeSet<Gene>();
		
		Iterator<? extends Annotation> iter = null;
		TreeSet<Annotation> exons = new TreeSet<Annotation>();
		exons.addAll(gene.getExonSet());
		boolean plusStrand = orientation.equals(Strand.POSITIVE);
		if(plusStrand) {
			iter = exons.iterator();
		} else {
			iter = exons.descendingIterator();
		}
		
		// Don't use first exon
		@SuppressWarnings("unused")
		Annotation firstExon = iter.next();
		int exonNumber = 1;
		
		while(iter.hasNext()) {
			Annotation exon = iter.next();
			if(exon.getSize() < numBases) {
				exonNumber++;
				continue;
			}
			int firstExonPos = plusStrand ? exon.getStart() : exon.getEnd() - 1;
			int windowFirstPos = plusStrand ? exon.getStart() : firstExonPos - numBases + 1;
			int windowLastPos = windowFirstPos + numBases - 1;
			int windowStart = Math.min(windowFirstPos, windowLastPos);
			int windowEnd = Math.max(windowFirstPos, windowLastPos) + 1;
			Gene window = new Gene(new BasicAnnotation(gene.getChr(), windowStart, windowEnd, orientation));
			window.setName(gene.getName() + "_downstream_of_splice_junction_" + exonNumber);
			rtrn.add(window);
			exonNumber++;
		}

		return rtrn;
		
	}

	@Override
	public boolean reverseDataIfMinusOrientation() {
		return true;
	}
	
}
