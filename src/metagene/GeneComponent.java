package metagene;

import java.util.Collection;
import java.util.Map;

import nextgen.core.annotation.Gene;

/**
 * @author prussell
 *
 */
public interface GeneComponent {

	/**
	 * Get name of gene component
	 * @return Name of gene component
	 */
	public String getName();
	
	/**
	 * Get the component
	 * @param gene The parent gene
	 * @return The component as a new annotation or null if not applicable
	 */
	public Collection<Gene> getComponent(Gene gene);

	/**
	 * Get the component for a set of genes
	 * @param genes Genes by chromosome
	 * @return The component by chromosome
	 */
	public Map<String, Collection<Gene>> getComponent(Map<String, Collection<Gene>> genes);
	
	/**
	 * Whether to reverse the order of data vector if region is on minus strand
	 * @return True if data should be reversed for negative strand annotations
	 */
	public boolean reverseDataIfMinusOrientation();
	
}
