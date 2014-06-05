package metagene;

import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import nextgen.core.annotation.Gene;


/**
 * @author prussell
 *
 */
public abstract class AbstractGeneComponent implements GeneComponent, Comparable<GeneComponent> {
	
	@Override
	public Map<String, Collection<Gene>> getComponent(Map<String, Collection<Gene>> genes) {
		Map<String, Collection<Gene>> rtrn = new TreeMap<String, Collection<Gene>>();
		for(String chr : genes.keySet()) {
			Collection<Gene> components = new TreeSet<Gene>();
			for(Gene gene : genes.get(chr)) {
				Collection<Gene> component = getComponent(gene);
				if(component != null) {
					components.addAll(component);
				}
			}
			if(!components.isEmpty()) {
				rtrn.put(chr, components);
			}
		}
		return rtrn;
	}
	
	@Override
	public int compareTo(GeneComponent other) {
		return getName().compareTo(other.getName());
	}
	
	@Override
	public boolean equals(Object o) {
		GeneComponent other = (GeneComponent)o;
		return getName().equals(other.getName());
	}
	
	@Override
	public int hashCode() {
		return getName().hashCode();
	}
	
}
