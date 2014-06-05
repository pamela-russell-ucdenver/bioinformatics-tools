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
public class GeneComponentIntron3PrimeEnd extends AbstractGeneComponent {

	private int length;
	
	/**
	 * @param numBases Number of bases to include
	 */
	public GeneComponentIntron3PrimeEnd(int numBases) {
		length = numBases;
	}
	
	@Override
	public String getName() {
		return "intron_3_prime_end_" + length;
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
		boolean plusStrand = orientation.equals(Strand.POSITIVE);
		
		Collection<Gene> rtrn = new TreeSet<Gene>();
		
		TreeSet<Annotation> introns = new TreeSet<Annotation>();
		introns.addAll(gene.getIntronSet());
		Iterator<? extends Annotation> iter = introns.iterator();
				
		int intronNumber = 0;
		
		while(iter.hasNext()) {
			Annotation intron = iter.next();
			if(intron.getSize() < length) {
				intronNumber++;
				continue;
			}
			int lastIntronPos = plusStrand ? intron.getEnd() - 1 : intron.getStart();
			int windowLastPos = plusStrand ? lastIntronPos : lastIntronPos + length - 1;
			int windowFirstPos = windowLastPos - length + 1;
			int windowStart = Math.min(windowFirstPos, windowLastPos);
			int windowEnd = Math.max(windowFirstPos, windowLastPos) + 1;
			Gene window = new Gene(new BasicAnnotation(gene.getChr(), windowStart, windowEnd, orientation));
			window.setName(gene.getName() + "_intron_3_prime_end_" + intronNumber);
			rtrn.add(window);
			intronNumber++;
		}

		return rtrn;
		
	}

	@Override
	public boolean reverseDataIfMinusOrientation() {
		return true;
	}
	
}
