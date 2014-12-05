package translation;

import net.sf.samtools.util.CloseableIterator;
import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotationcollection.AnnotationCollection;
import guttmanlab.core.annotationcollection.FeatureCollection;

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
	
	/**
	 * Get all upstream ORFs of a gene
	 * @param orfFinder ORF finder object
	 * @param gene The gene
	 * @return Collection of uORFs
	 */
	public static AnnotationCollection<UpstreamORF> findAllUpstreamORFs(ORFFinder orfFinder, Gene gene) {
		Annotation utr5 = gene.get5UTR();
		FeatureCollection<UpstreamORF> rtrn = new FeatureCollection<UpstreamORF>(orfFinder.getCoordSpace());
		if(utr5 == null) {
			return rtrn;
		}
		AnnotationCollection<Gene> uorfs = orfFinder.getAllORFs(utr5);
		CloseableIterator<Gene> iter = uorfs.sortedIterator();
		while(iter.hasNext()) {
			Gene uorf = iter.next();
			int cdsStart = uorf.getCodingRegion().getReferenceStartPosition();
			int cdsEnd = uorf.getCodingRegion().getReferenceEndPosition();
			rtrn.addAnnotation(new UpstreamORF(new Gene(gene.getBlockSet(), cdsStart, cdsEnd, uorf.getName()), gene));
		}
		return rtrn;
	}
	
}
