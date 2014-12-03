package translation;

import java.util.Collection;

import guttmanlab.core.annotation.Gene;

public class UpstreamORF extends Gene {

	private Gene parentGene;
	
	/**
	 * @param orf The uORF as a gene object
	 * @param parent The parent gene
	 */
	public UpstreamORF(Gene orf, Gene parent) {
		super(orf.getBlockSet(), orf.getCodingRegion().getReferenceStartPosition(), orf.getCodingRegion().getReferenceEndPosition(), orf.getName());
		parentGene = parent;
	}

	/**
	 * @return Parent gene
	 */
	public Gene getParent() {
		return parentGene;
	}
	
	public static Collection<UpstreamORF> findAllUpstreamORFs(Gene gene) {
		
		
		
	}
	
}
