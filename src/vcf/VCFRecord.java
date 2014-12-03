package vcf;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.sequence.Sequence;
import guttmanlab.core.util.StringParser;

import java.util.ArrayList;
import java.util.List;

/**
 * A VCF record
 * @author prussell
 *
 */
public class VCFRecord {
	
	/**
	 * The alternate allele field
	 * @author prussell
	 *
	 */
	private class AlternateAllele {
		
		private List<Sequence> alternateSeqs;
		
		/**
		 * @param vcfField The string alternate allele field from VCF line
		 */
		public AlternateAllele(String vcfField) {
			
			if(vcfField.contains("[") || vcfField.contains("]")) {
				throw new UnsupportedOperationException("Breakends not supported");
			}
			
			if(vcfField.contains("<") || vcfField.contains(">")) {
				throw new UnsupportedOperationException("ID strings not supported");
			}
			
			StringParser s = new StringParser();
			s.parse(vcfField, ",");
			alternateSeqs = new ArrayList<Sequence>();
			for(int i = 0; i < s.getFieldCount(); i++) {
				alternateSeqs.add(new Sequence(s.asString(i)));
			}
			
		}
		
		public String toString() {
			String rtrn = alternateSeqs.get(0).getSequenceBases();
			for(int i = 1; i < alternateSeqs.size(); i++) {
				rtrn += "," + alternateSeqs.get(i).getSequenceBases();
			}
			return rtrn;
		}
		
		/**
		 * Reverse complement the alleles
		 */
		public void reverseComplement() {
			List<Sequence> rc = new ArrayList<Sequence>();
			for(int i = 0; i < alternateSeqs.size(); i++) {
				rc.add(alternateSeqs.get(i).reverseComplement());
			}
			alternateSeqs = rc;
		}
		
	}
	
	private String chrom;
	private int pos; // One based
	private String id;
	private Sequence ref;
	private AlternateAllele alt;
	private double qual;
	private String filter;
	private String info;
	private List<String> genotypes;
	
	private VCFRecord() {}
	
	/**
	 * Make a copy
	 * @param other Other record
	 * @return Copy of other record
	 */
	public static VCFRecord copy(VCFRecord other) {
		VCFRecord rtrn = new VCFRecord();
		rtrn.chrom = other.chrom;
		rtrn.pos = other.pos;
		rtrn.id = other.id;
		rtrn.ref = other.ref;
		rtrn.alt = other.alt;
		rtrn.qual = other.qual;
		rtrn.filter = other.filter;
		rtrn.info = other.info;
		rtrn.genotypes = other.genotypes;
		return rtrn;
	}
	
	/**
	 * @param line VCF file line
	 */
	public VCFRecord(String line) {
		StringParser s = new StringParser();
		s.parse(line);
		chrom = s.asString(0);
		pos = s.asInt(1);
		id = s.asString(2);
		ref = new Sequence(s.asString(3));
		alt = new AlternateAllele(s.asString(4));
		qual = s.asDouble(5);
		filter = s.asString(6);
		info = s.asString(7);
		genotypes = new ArrayList<String>();
		for(int i = 8; i < s.getFieldCount(); i++) {
			genotypes.add(s.asString(i));
		}
	}
	
	/**
	 * Get the position as zero based
	 * @return Zero based position
	 */
	public int getZeroBasedPos() {
		return pos - 1;
	}
	
	/**
	 * Get the position as one based
	 * @return One based position
	 */
	public int getOneBasedPos() {
		return pos;
	}
	
	/**
	 * Set position
	 * @param zeroBasedPos Zero based position to set
	 */
	private void setZeroBasedPos(int zeroBasedPos) {
		pos = zeroBasedPos + 1;
	}
	
	/**
	 * Set position
	 * @param oneBasedPos One based position to set
	 */
	@SuppressWarnings("unused")
	private void setOneBasedPos(int oneBasedPos) {
		pos = oneBasedPos;
	}
	
	/**
	 * Reverse complement the reference and alternate alleles
	 */
	private void reverseComplementAlleles() {
		ref = ref.reverseComplement();
		alt.reverseComplement();
	}
	
	/**
	 * Set chromosome
	 * @param newChrom New chromosome
	 */
	private void setChrom(String newChrom) {
		chrom = newChrom;
	}
	
	public String toString() {
		String rtrn = chrom + "\t" + pos + "\t" + id + "\t" + ref.getSequenceBases() + "\t" + alt.toString() + "\t" + qual + "\t" + filter + "\t" + info;
		for(String genotype : genotypes) {
			rtrn += "\t" + genotype;
		}
		return rtrn;
	}
	
	public String getChrom() {
		return chrom;
	}
	
	/**
	 * Convert the entire record to transcript coordinates
	 * @param record The record
	 * @param transcript The transcript
	 * @return The converted record or null if position doesn't overlap transcript
	 */
	public static VCFRecord convertToTranscriptCoords(VCFRecord record, Annotation transcript) {
		VCFRecord rtrn = copy(record);
		Strand strand = transcript.getOrientation();
		if(strand.equals(Strand.BOTH) || strand.equals(Strand.INVALID) || strand.equals(Strand.UNKNOWN)) {
			throw new IllegalArgumentException("Annotation strand must be known.");
		}
		// Reverse complement alleles if necessary
		if(strand.equals(Strand.NEGATIVE)) {
			rtrn.reverseComplementAlleles();
		}
		// Set "chromosome" to transcript name
		rtrn.setChrom(transcript.getName());
		// Set pos to transcript coordinate
		int zeroBasedGenomeCoord = rtrn.getZeroBasedPos();
		int zeroBasedTranscriptCoord = transcript.getRelativePositionFrom5PrimeOfFeature(zeroBasedGenomeCoord);
		if(zeroBasedTranscriptCoord < 0) { // Coordinate doesn't overlap transcript
			return null;
		}
		rtrn.setZeroBasedPos(zeroBasedTranscriptCoord);
		return rtrn;
	}
	
}
