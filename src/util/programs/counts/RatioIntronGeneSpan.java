package util.programs.counts;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.Optional;

import org.apache.log4j.Logger;

import net.sf.samtools.util.CloseableIterator;
import guttmanlab.core.annotation.BEDFileRecord;
import guttmanlab.core.annotation.BlockedAnnotation;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.AnnotationCollection;
import guttmanlab.core.coordinatespace.CoordinateSpace;
import guttmanlab.core.util.CommandLineParser;
import guttmanlab.core.util.CountLogger;
import htsjdk.samtools.fork.CigarElement;
import htsjdk.samtools.fork.CigarOperator;
import htsjdk.samtools.fork.SAMRecord;
import htsjdk.samtools.fork.SAMRecordIterator;
import htsjdk.samtools.fork.SamReader;
import htsjdk.samtools.fork.SamReaderFactory;

/**
 * Measure proportion of reads mapping to introns and exons
 * @author prussell
 *
 */
public class RatioIntronGeneSpan {
	
	
	private static class TranscriptStats {
		
		private String id;
		
		public TranscriptStats(String transcriptID) {
			id = transcriptID;
		}
		
		@SuppressWarnings("unused")
		public int getNumExons() {
			return numExons.get().intValue();
		}
		
		public void setNumExons(int numExons) {
			this.numExons = Optional.of(Integer.valueOf(numExons));
		}
		
		@SuppressWarnings("unused")
		public int getTotalExonSize() {
			return totalExonSize.get().intValue();
		}
		
		public void setTotalExonSize(int totalExonSize) {
			this.totalExonSize = Optional.of(Integer.valueOf(totalExonSize));
		}
		
		@SuppressWarnings("unused")
		public int getGeneSpan() {
			return geneSpan.get().intValue();
		}
		
		public void setGeneSpan(int geneSpan) {
			this.geneSpan = Optional.of(Integer.valueOf(geneSpan));
		}
		
		public long getIntronCount() {
			return intronCount.get().longValue();
		}
		
		/**
		 * Get count of reads not contained in introns PLUS
		 * reads contained in introns that are SPLICED
		 * @return
		 */
		public Optional<Long> getCountExonPlusSpliced() {
			try {
				return Optional.of(Long.valueOf(getGeneSpanCount() - getIntronCount()));
			} catch(NoSuchElementException e) {
				return Optional.empty();
			}
		}
		
		public long getGeneSpanCount() {
			return geneSpanCount.get().longValue();
		}
		
		public void setIntronCount(long intronCount) {
			this.intronCount = Optional.of(Long.valueOf(intronCount));
		}
		
		public void setGeneSpanCount(long geneSpanCount) {
			this.geneSpanCount = Optional.of(Long.valueOf(geneSpanCount));
		}
		
		
		public Optional<Float> getProportionOfReadsFromExons() {
			Optional<Float> introns = getProportionOfReadsFromIntrons();
			if(introns.isPresent()) return Optional.of(Float.valueOf(1 - introns.get().floatValue()));
			return Optional.empty();
		}
		
		public Optional<Float> getProportionOfReadsFromIntrons() {
			if(intronCount.isPresent() && geneSpanCount.isPresent()) {
				long g = geneSpanCount.get();
				long e = intronCount.get();
				if(g != 0) return Optional.of(Float.valueOf((float)e / g));
				else return Optional.empty();
			}
			throw new NoSuchElementException();
		}
		
		public static final String columnHeaders() {
			return 	"transcript_id\t"
					+"gene_span_length\t"
					+ "num_exons\t"
					+ "total_length_of_exons\t"
					+ "total_count_over_gene_span\t"
					+ "total_count_over_exons_plus_all_spliced\t"
					+ "total_count_over_introns\t"
					+ "proportion_exons\t"
					+ "proportion_introns\t";
		}
		
		public String toString() {
			String rtrn = id + "\t";
			if(geneSpan.isPresent()) rtrn += geneSpan.get() + "\t";
			else rtrn += "-\t";
			if(numExons.isPresent()) rtrn += numExons.get() + "\t";
			else rtrn += "-\t";
			if(totalExonSize.isPresent()) rtrn += totalExonSize.get() + "\t";
			else rtrn += "-\t";
			if(geneSpanCount.isPresent()) rtrn += geneSpanCount.get() + "\t";
			else rtrn += "-\t";
			if(getCountExonPlusSpliced().isPresent()) rtrn += getCountExonPlusSpliced().get() + "\t";
			else rtrn += "-\t";
			if(intronCount.isPresent()) rtrn += intronCount.get() + "\t";
			else rtrn += "-\t";
			try {
				rtrn += getProportionOfReadsFromExons().get() + "\t";
				rtrn += getProportionOfReadsFromIntrons().get() + "\t";
			} catch(NoSuchElementException e) {
				rtrn += "-\t-\t";
			}
			return rtrn;
		}
		
		private Optional<Integer> numExons;
		private Optional<Integer> totalExonSize;
		private Optional<Integer> geneSpan;
		private Optional<Long> intronCount;
		private Optional<Long> geneSpanCount;
		
	}
	
	private static Logger logger = Logger.getLogger(RatioIntronGeneSpan.class.getName());
	
	/**
	 * Get count of UNSPLICED reads that are FULLY CONTAINED in introns
	 * @param transcript Transcript
	 * @param samReader SAM reader
	 * @return Count of UNSPLICED reads that are FULLY CONTAINED in introns.
	 * Spliced reads that are nevertheless contained in the intron are NOT COUNTED.
	 */
	private static long getTotalCountOverIntrons(BlockedAnnotation transcript, SamReader samReader) {
		String chr = transcript.getReferenceName();
		Collection<SingleInterval> introns = new ArrayList<SingleInterval>();
		Iterator<SingleInterval> blockIter = transcript.getBlocks();
		int lastBlockEnd = blockIter.next().getReferenceEndPosition();
		while(blockIter.hasNext()) {
			SingleInterval nextBlock = blockIter.next();
			int nextBlockStart = nextBlock.getReferenceStartPosition();
			introns.add(new SingleInterval(chr, lastBlockEnd, nextBlockStart));
			lastBlockEnd = nextBlock.getReferenceEndPosition();
		}
		return introns
				.stream()
				.mapToLong(block -> getCount(chr, block.getReferenceStartPosition(), block.getReferenceEndPosition(), samReader, true))
				.sum();
	}
	
	private static long getCount(String chr, int start, int end, SamReader samReader, boolean unsplicedOnly) {
		SAMRecordIterator iter = samReader.queryContained(chr, start, end);
		long rtrn = 0;
		while(iter.hasNext()) {
			SAMRecord record = iter.next();
			if(unsplicedOnly) {
				Iterator<CigarElement> cigarIter = record.getCigar().iterator();
				int numBlocks = 0;
				while(cigarIter.hasNext()) {
					if(cigarIter.next().getOperator().equals(CigarOperator.M)) {numBlocks++;}
				}
				if(numBlocks > 1) continue;
			}
			if(!record.getReadUnmappedFlag()) rtrn++;
		}
		iter.close();
		//if(rtrn > 0) logger.info("COUNT\tunspliced_only:" + unsplicedOnly + "\t" + chr + ":" + start + "-" + end + "\t" + rtrn);
		return rtrn;
	}
	
	private static long getCountOverGeneSpan(BlockedAnnotation transcript, SamReader samReader) {
		return getCount(transcript.getReferenceName(), transcript.getReferenceStartPosition(), 
				transcript.getReferenceEndPosition(), samReader, false);
	}
	
	private static int getTotalSizeOfExons(BlockedAnnotation transcript) {
		return transcript.getBlockSet().stream()
				.mapToInt(block -> block.size())
				.sum();
	}
	
	private static TranscriptStats getTranscriptStats(BlockedAnnotation transcript, SamReader samReader) {
		TranscriptStats rtrn = new TranscriptStats(transcript.getName());
		rtrn.setIntronCount(getTotalCountOverIntrons(transcript, samReader));
		rtrn.setGeneSpan(transcript.getReferenceEndPosition() - transcript.getReferenceStartPosition());
		rtrn.setGeneSpanCount(getCountOverGeneSpan(transcript, samReader));
		rtrn.setNumExons(transcript.getNumberOfBlocks());
		rtrn.setTotalExonSize(getTotalSizeOfExons(transcript));
		return rtrn;
	}
	
	private static void writeTranscriptStats(BlockedAnnotation transcript, SamReader samReader, FileWriter writer) {
		try {
			TranscriptStats stats = getTranscriptStats(transcript, samReader);
			if(stats.getGeneSpanCount() == 0) return;
			writer.write(stats.toString() + "\n");
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	private static void writeTranscriptStats(AnnotationCollection<? extends BlockedAnnotation> transcripts,
			SamReader samReader, FileWriter writer) {
		CountLogger countLogger = new CountLogger(transcripts.getNumAnnotations(), 10);
		CloseableIterator<? extends BlockedAnnotation> iter = transcripts.sortedIterator();
		while(iter.hasNext()) {
			countLogger.advance();
			writeTranscriptStats(iter.next(), samReader, writer);
		}
		iter.close();
	}
	
	private static void writeTranscriptStats(File bedFile, CoordinateSpace coordSpace, File bamFile, File outputTable) {
		try {
			AnnotationCollection<BEDFileRecord> transcripts = BEDFileIO.loadFromFile(bedFile, coordSpace);
			FileWriter writer = new FileWriter(outputTable);
			SamReader samReader = SamReaderFactory.makeDefault().open(bamFile);
			writer.write(TranscriptStats.columnHeaders() + "\n");
			writeTranscriptStats(transcripts, samReader, writer);
			writer.close();
			samReader.close();
		} catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	public static void main(String[] args) {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-b", "Bam file", true);
		p.addStringArg("-t", "Bed file", true);
		p.addStringArg("-c", "Chromosome size file", true);
		p.addStringArg("-o", "Output table", true);
		p.parse(args);
		File bedFile = new File(p.getStringArg("-t"));
		CoordinateSpace coordSpace = new CoordinateSpace(p.getStringArg("-c"));
		File bamFile = new File(p.getStringArg("-b"));
		File output = new File(p.getStringArg("-o"));
		
		writeTranscriptStats(bedFile, coordSpace, bamFile, output);
		
		logger.info("");
		logger.info("All done");
		
	}
	
}
