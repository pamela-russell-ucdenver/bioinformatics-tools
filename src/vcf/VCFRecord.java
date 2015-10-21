package vcf;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.sequence.Sequence;
import guttmanlab.core.util.StringParser;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.apache.log4j.Logger;

/**
 * A VCF record
 * @author prussell
 *
 */
public class VCFRecord {
	
	public static Logger logger = Logger.getLogger(VCFRecord.class.getName());
	
	/**
	 * The alternate allele field
	 * @author prussell
	 *
	 */
	public class AlternateAllele {
		
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
		
		/**
		 * @param singleAltSequence A single alternate sequence
		 */
		public AlternateAllele(Sequence singleAltSequence) {
			alternateSeqs = new ArrayList<Sequence>();
			alternateSeqs.add(singleAltSequence);
		}
		
		/**
		 * If there is one alternate sequence, get it
		 * Otherwise throw an exception
		 * @return The only alternate sequence
		 */
		public Sequence getSingleSeq() {
			if(alternateSeqs.size() != 1) {
				throw new IllegalStateException("There is not exactly one alternate sequence");
			}
			return alternateSeqs.iterator().next();
		}
		
		/**
		 * If there is one alternate sequence, get its size
		 * Otherwise throw an exception
		 * @return The size of the only alternate sequence
		 */
		public int getSizeOfSingleSeq() {
			if(alternateSeqs.size() != 1) {
				throw new IllegalStateException("There is not exactly one alternate sequence");
			}
			return alternateSeqs.iterator().next().getSequenceBases().length();
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
		public AlternateAllele reverseComplement() {
			String list = "";
			for(int i = 0; i < alternateSeqs.size(); i++) {
				if(i > 0) {
					list += ",";
				}
				list += alternateSeqs.get(i).reverseComplement().getSequenceBases();
			}
			return new AlternateAllele(list);
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
	private String format;
	private List<String> sampleIDs;
	private Map<String, String> genotypes;
	
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
	 * @param headerLine Header line beginning with #CHROM
	 */
	public VCFRecord(String line, String headerLine) {
		loadHeader(headerLine);
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
		genotypes = new HashMap<String, String>();
		if(sampleIDs != null) {
			format = s.asString(8);
			if(s.getFieldCount() != sampleIDs.size() + 9) {
				throw new IllegalArgumentException("Line must have 9 more fields (" + s.getFieldCount() +  ") than sample IDs (" + sampleIDs.size() + ")");
			}
			for(int i = 0; i < sampleIDs.size(); i++) {
				genotypes.put(sampleIDs.get(i), s.asString(i + 9));
			}
		}
	}
	
	private void loadHeader(String headerLine) {
		StringParser p = new StringParser();
		p.parse(headerLine);
		if(p.getFieldCount() < 9) return;
		sampleIDs = new ArrayList<String>();
		for(int i = 9; i < p.getFieldCount(); i++) {
			sampleIDs.add(p.asString(i));
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
	 * Get the reference allele as a sequence object
	 * @return
	 */
	public Sequence getRefAllele() {
		return ref;
	}
	
	/**
	 * Get the alternate allele
	 * @return The alternate allele
	 */
	public AlternateAllele getAlternateAllele() {
		return alt;
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
		alt = alt.reverseComplement();
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
		if(sampleIDs != null) {
			rtrn += "\t" + format;
			for(String sample : sampleIDs) {
				rtrn += "\t" + genotypes.get(sample);
			}
		}
		return rtrn;
	}
	
	/**
	 * @return Chromosome name
	 */
	public String getChrom() {
		return chrom;
	}
	
	/**
	 * Get a single interval representing the reference allele region
	 * Orientation is set to BOTH
	 * @return Annotation representing the reference allele
	 */
	public Annotation getRefAlleleAsAnnotation() {
		int end = pos + ref.getSequenceBases().length();
		return new SingleInterval(chrom, pos, end, Strand.BOTH);
	}
	
	/**
	 * Convert the entire record to transcript coordinates
	 * @param record The record
	 * @param transcript The transcript
	 * @return The converted record or null if position doesn't overlap transcript
	 */
	public static VCFRecord convertToTranscriptCoords(VCFRecord record, Annotation transcript, Sequence transcriptSequence) {
		// Make sure ref allele is fully contained within a single block of the transcript
		Annotation refAllele = record.getRefAlleleAsAnnotation();
		int refStart = refAllele.getReferenceStartPosition();
		int refEnd = refAllele.getReferenceEndPosition();
		boolean contains = false;
		Iterator<SingleInterval> blockIter = transcript.getBlocks();
		while(blockIter.hasNext()) {
			SingleInterval block = blockIter.next();
			if(refStart >= block.getReferenceStartPosition() && refEnd <= block.getReferenceEndPosition()) {
				contains = true;
				break;
			}
		}
		if(!contains) {
			//logger.info(transcript.getName() + " does not contain ref allele " + record.getRefAlleleAsAnnotation().toUCSC());
			return null;
		}
		VCFRecord rtrn = copy(record);
		Strand strand = transcript.getOrientation();
		if(strand.equals(Strand.BOTH) || strand.equals(Strand.INVALID) || strand.equals(Strand.UNKNOWN)) {
			throw new IllegalArgumentException("Annotation strand must be known.");
		}
		// Set "chromosome" to transcript name
		rtrn.setChrom(transcript.getName());
		// Set pos to transcript coordinate
		int zeroBasedGenomeCoord = rtrn.getZeroBasedPos();
		int zeroBasedTranscriptCoord = transcript.getRelativePositionFrom5PrimeOfFeature(zeroBasedGenomeCoord);
		
		/*
		 *  Special procedure for deletion on negative strand gene
		 *  Identify region to be deleted
		 *  "Ref" allele is one upstream base plus deletion region
		 *  "Alt" allele is just the one upstream base
		 *  Get upstream base from transcript sequence
		 */
		int refSize = refAllele.size();
		int altSize =  record.getAlternateAllele().getSizeOfSingleSeq();
		if(refSize > 1 && altSize == 1 && strand.equals(Strand.NEGATIVE)) {
			// Start of actual deleted region on genome is zeroBasedGenomeCoord + 1
			// Deletion size is refSize - 1
			zeroBasedTranscriptCoord = transcript.getRelativePositionFrom5PrimeOfFeature(zeroBasedGenomeCoord + 1) - refSize + 1;
			String upstreamBaseToKeep = transcriptSequence.getSequenceBases().substring(zeroBasedTranscriptCoord, zeroBasedTranscriptCoord + 1).toUpperCase();
			Sequence deletedRegionGenomic = new Sequence(record.getRefAllele().getSequenceBases().substring(1));
			String deletedRegionTranscript = deletedRegionGenomic.reverseComplement().getSequenceBases().toUpperCase();
			String transcriptRef = upstreamBaseToKeep + deletedRegionTranscript;
			// Set fields in return record
			rtrn.setZeroBasedPos(zeroBasedTranscriptCoord);
			rtrn.ref = new Sequence(transcriptRef);
			rtrn.alt = new VCFRecord().new AlternateAllele(new Sequence(upstreamBaseToKeep));			
			return rtrn;
		}
		
		/*
		 *  Special procedure for insertion on negative strand gene
		 *  Identify region to be inserted
		 *  "Alt" allele is one upstream base plus insertion region
		 *  "Ref" allele is just the one upstream base
		 *  Get upstream base from transcript sequence
		 */
		if(refSize == 1 && altSize > 1 && strand.equals(Strand.NEGATIVE)) {
			zeroBasedTranscriptCoord = transcript.getRelativePositionFrom5PrimeOfFeature(zeroBasedGenomeCoord + 1);
			String upstreamBaseToKeep = transcriptSequence.getSequenceBases().substring(zeroBasedTranscriptCoord, zeroBasedTranscriptCoord + 1).toUpperCase();
			Sequence insertedSequenceGenomic = new Sequence(record.getAlternateAllele().getSingleSeq().getSequenceBases().substring(1));
			String insertedSequenceTranscript = insertedSequenceGenomic.reverseComplement().getSequenceBases().toUpperCase();
			String transcriptAlt = upstreamBaseToKeep + insertedSequenceTranscript;
			rtrn.setZeroBasedPos(zeroBasedTranscriptCoord);
			rtrn.alt = new VCFRecord().new AlternateAllele(new Sequence(transcriptAlt));
			rtrn.ref = new Sequence(upstreamBaseToKeep);
			return rtrn;
		}
		
		
		
		if(zeroBasedTranscriptCoord < 0) { // Coordinate doesn't overlap transcript
			return null;
		}
		// Reverse complement alleles if necessary
		if(strand.equals(Strand.NEGATIVE)) {
			rtrn.reverseComplementAlleles();
		}
		rtrn.setZeroBasedPos(zeroBasedTranscriptCoord);
		return rtrn;
	}
	
}
