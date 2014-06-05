package metagene;

import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.log4j.Logger;

import nextgen.core.annotation.Gene;

/**
 * @author prussell
 *
 */
public class GenePartition {
	
	private Collection<GeneComponent> components;
	private static Logger logger = Logger.getLogger(GenePartition.class.getName());
	
	/**
	 * @param geneComponents The set of gene components for the partition
	 */
	public GenePartition(Collection<GeneComponent> geneComponents) {
		if(geneComponents.isEmpty()) {
			throw new IllegalArgumentException("Set of gene components cannot be empty.");
		}
		components = geneComponents;
	}
	
	/**
	 * Get the set of components for the partition
	 * @return The set of components
	 */
	public Collection<GeneComponent> getComponents() {
		return components;
	}
	
	/**
	 * Get the partition of the gene
	 * @param gene The gene
	 * @return Map of component type to component of the gene
	 */
	public Map<GeneComponent, Collection<Gene>> getPartition(Gene gene) {
		Map<GeneComponent, Collection<Gene>> rtrn = new TreeMap<GeneComponent, Collection<Gene>>();
		for(GeneComponent component : components) {
			rtrn.put(component, component.getComponent(gene));
		}
		return rtrn;
	}
	
	/**
	 * Get the partition of all the genes together
	 * @param genes Genes by chromosome
	 * @return Map of component type to set of the component over all genes
	 */
	public Map<GeneComponent, Collection<Gene>> getPartition(Map<String, Collection<Gene>> genes) {
		logger.info("");
		logger.info("Getting full gene partition...");
		Map<GeneComponent, Collection<Gene>> rtrn = new TreeMap<GeneComponent, Collection<Gene>>();
		for(GeneComponent component : components) {
			Collection<Gene> thisComponent = new TreeSet<Gene>();
			for(String chr : genes.keySet()) {
				logger.info(chr + "\t" + component.getName());
				for(Gene gene : genes.get(chr)) {
					Collection<Gene> componentThisGene = component.getComponent(gene);
					if(componentThisGene == null) {
						continue;
					}
					for(Gene c : componentThisGene) {
						if(c != null) {
							thisComponent.add(c);
						}
					}
				}
			}
			rtrn.put(component, thisComponent);
		}
		return rtrn;		
	}
	
}
