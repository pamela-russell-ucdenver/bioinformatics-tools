package metagene;

import java.util.ArrayList;
import java.util.Collection;
import java.util.TreeSet;

import nextgen.core.annotation.Gene;

/**
 * @author prussell
 *
 */
public class GeneComponentCDS extends AbstractGeneComponent {
	
	private static String NAME = "CDS";
	private Collection<GeneComponent> requiredComponents;
	
	/**
	 * 
	 */
	public GeneComponentCDS() {
		this(null);
	}
	
	/**
	 * @param requiredGeneComponents Components that gene must have in order to report a CDS
	 */
	public GeneComponentCDS(Collection<GeneComponent> requiredGeneComponents) {
		requiredComponents = new ArrayList<GeneComponent>();
		if(requiredGeneComponents != null) {
			requiredComponents.addAll(requiredGeneComponents);
		}
	}

	@Override
	public String getName() {
		return NAME;
	}

	@Override
	public Collection<Gene> getComponent(Gene gene) {
		if(!gene.hasCDS()) {
			return null;
		}
		for(GeneComponent requirement : requiredComponents) {
			if(requirement.getComponent(gene) == null) {
				return null;
			}
		}
		Collection<Gene> rtrn = new TreeSet<Gene>();
		rtrn.add(gene.getCDS());
		return rtrn;
	}

	@Override
	public boolean reverseDataIfMinusOrientation() {
		return true;
	}

	
}
