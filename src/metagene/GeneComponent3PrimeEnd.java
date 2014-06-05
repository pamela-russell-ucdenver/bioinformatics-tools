/**
 * 
 */
package metagene;

import java.util.Collection;
import java.util.TreeSet;

import nextgen.core.annotation.Annotation.Strand;
import nextgen.core.annotation.Gene;

/**
 * @author prussell
 *
 */
public class GeneComponent3PrimeEnd extends AbstractGeneComponent {
	
	private int size;
	
	/**
	 * @param length Window size
	 */
	public GeneComponent3PrimeEnd(int length) {
		size = length;
	}

	/* (non-Javadoc)
	 * @see metagene.GeneComponent#getName()
	 */
	@Override
	public String getName() {
		return "3_prime_end";
	}

	/* (non-Javadoc)
	 * @see metagene.GeneComponent#getComponent(nextgen.core.annotation.Gene)
	 */
	@Override
	public Collection<Gene> getComponent(Gene gene) {
		Strand strand = gene.getOrientation();
		if(strand.equals(Strand.UNKNOWN)) {
			return null;
		}
		if(gene.size() < size) {
			return null;
		}
		Collection<Gene> rtrn = new TreeSet<Gene>();
		
		int windowLastPos = gene.transcriptToGenomicPosition(gene.size() - 1);
		int windowFirstPos = gene.transcriptToGenomicPosition(gene.size() - size);
		int windowStart = Math.min(windowFirstPos, windowLastPos);
		int windowEnd = Math.max(windowFirstPos, windowLastPos) + 1;
		Gene window = gene.copy();
		window.trimAbsolute(windowStart, windowEnd);
		window.setName(gene.getName() + "_3_prime_end");
		rtrn.add(window);

		return rtrn;

	}

	/* (non-Javadoc)
	 * @see metagene.GeneComponent#reverseDataIfMinusOrientation()
	 */
	@Override
	public boolean reverseDataIfMinusOrientation() {
		return true;
	}

}
