package metagene;

import java.util.Collection;
import java.util.TreeSet;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Gene;

/**
 * @author prussell
 *
 */
public class GeneComponentIntrons extends AbstractGeneComponent {

	@Override
	public String getName() {
		return "introns";
	}

	@Override
	public Collection<Gene> getComponent(Gene gene) {
		Collection<? extends Annotation> introns = gene.getIntronSet();
		Collection<Gene> rtrn = new TreeSet<Gene>();
		int intronNumber = 1;
		for(Annotation intron : introns) {
			Gene g = new Gene(intron);
			g.setName(gene.getName() + "_intron_" + intronNumber);
			intronNumber++;
			rtrn.add(g);
		}
		return rtrn;
	}

	@Override
	public boolean reverseDataIfMinusOrientation() {
		return true;
	}

	
}
