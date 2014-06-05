package metagene;

import java.util.Collection;
import java.util.TreeSet;

import nextgen.core.annotation.Gene;

/**
 * @author prussell
 *
 */
public class GeneComponentWholeGene extends AbstractGeneComponent {

	@Override
	public String getName() {
		return "entire_gene";
	}

	@Override
	public Collection<Gene> getComponent(Gene gene) {
		Collection<Gene> rtrn = new TreeSet<Gene>();
		rtrn.add(gene);
		return rtrn;
	}
	
	@Override
	public boolean reverseDataIfMinusOrientation() {
		return true;
	}


}
