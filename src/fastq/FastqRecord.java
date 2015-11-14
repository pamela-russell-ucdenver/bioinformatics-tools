package fastq;

import java.io.BufferedWriter;
import java.io.IOException;

import broad.core.sequence.Sequence;
import net.sf.samtools.SAMRecord;

public class FastqRecord {
	
	private String sequence;
	private String quality;
	private String name;
	private String description;
	
	/**
	 * @param seq Read sequence
	 * @param qual Quality string
	 * @param id Read name
	 */
	public FastqRecord(String seq, String qual, String id) {
		this(seq, qual, id, "");
	}
	
	/**
	 * @param seq Read sequence
	 * @param qual Quality string
	 * @param id Read name
	 * @param desc Description
	 */
	public FastqRecord(String seq, String qual, String id, String desc) {
		sequence = seq;
		quality = qual;
		name = id;
		description = desc;
	}
	
	/**
	 * Get a fastq record by reverting the SAM record
	 * If record is mapped on minus strand, reverse complement SEQ and reverse QUAL
	 * No description
	 * @param samRecord SAM record
	 */
	public FastqRecord(SAMRecord samRecord) {
		this(revertSamSequence(samRecord), revertSamQual(samRecord), samRecord.getReadName());
	}
	
	/**
	 * Get a fastq record by reverting the SAM record
	 * If record is mapped on minus strand, reverse complement SEQ and reverse QUAL
	 * No description
	 * @param samRecord SAM record
	 * @param name Name to use instead of QNAME from SAM record
	 */
	public FastqRecord(SAMRecord samRecord, String name) {
		this(revertSamSequence(samRecord), revertSamQual(samRecord), name);
	}
	
	private static String revertSamSequence(SAMRecord samRecord) {
		String alignedStrandSeq = samRecord.getReadString();
		if(samRecord.getReadNegativeStrandFlag()) {
			return Sequence.reverseSequence(alignedStrandSeq);
		} else {
			return alignedStrandSeq;
		}
	}
	
	private static String revertSamQual(SAMRecord samRecord) {
		String alignedStrandQual = samRecord.getBaseQualityString();
		if(samRecord.getReadNegativeStrandFlag()) {
			return new StringBuilder(alignedStrandQual).reverse().toString();
		} else {
			return alignedStrandQual;
		}
	}
	
	public void write(BufferedWriter w) throws IOException {
		w.write("@" + name);
		w.newLine();
		w.write(sequence);
		w.newLine();
		w.write("+" + description);
		w.newLine();
		w.write(quality);
		w.newLine();
	}

}
