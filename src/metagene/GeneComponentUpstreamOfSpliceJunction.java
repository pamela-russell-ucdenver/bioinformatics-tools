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
public class GeneComponentUpstreamOfSpliceJunction extends AbstractGeneComponent {

	private int numBases;
	
	/**
	 * @param numBasesUpstream Number of bases to include
	 */
	public GeneComponentUpstreamOfSpliceJunction(int numBasesUpstream) {
		numBases = numBasesUpstream;
	}
	
	@Override
	public String getName() {
		return "upstream_of_splice_junction_" + numBases;
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
		int exonNumber = 1;
		while(iter.hasNext()) {
			Annotation exon = iter.next();
			if(!iter.hasNext()) {
				// Don't use the last exon
				continue;
			}
			if(exon.getSize() < numBases) {
				exonNumber++;
				continue;
			}
			int lastExonPos = plusStrand ? exon.getEnd() - 1 : exon.getStart();
			int windowLastPos = plusStrand ? lastExonPos : lastExonPos + numBases - 1;
			int windowFirstPos = windowLastPos - numBases + 1;
			int windowStart = Math.min(windowFirstPos, windowLastPos);
			int windowEnd = Math.max(windowFirstPos, windowLastPos) + 1;
			Gene window = new Gene(new BasicAnnotation(gene.getChr(), windowStart, windowEnd, orientation));
			window.setName(gene.getName() + "_upstream_of_splice_junction_" + exonNumber);
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
