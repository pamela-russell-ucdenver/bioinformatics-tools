package util.programs.counts;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Optional;
import java.util.function.Predicate;

import org.apache.log4j.Logger;

import net.sf.samtools.util.CloseableIterator;
import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Annotation.Strand;
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
import htsjdk.samtools.fork.AlignmentBlock;

/**
 * Measure proportion of reads mapping to introns and exons
 * @author prussell
 *
 */
public class RatioIntronGeneSpan {
	
	private AnnotationCollection<BEDFileRecord> transcriptome;
	private SamReader samReader;
	private String sampleId;
	
	private RatioIntronGeneSpan(File bamFile, File transcriptomeBed, CoordinateSpace coordSpace, String sampleId) {
		this.sampleId = sampleId;
		try {
			transcriptome = BEDFileIO.loadFromFile(transcriptomeBed, coordSpace);
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		samReader = SamReaderFactory.makeDefault().open(bamFile);
	}
	
	
	
	private class TranscriptStats {
		
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
		
		@SuppressWarnings("unused")
		public long getIntronCountUnspliced() {
			return intronCountUnspliced.get().longValue();
		}
		
		public long getIntronCountAll() {
			return intronCountAll.get().longValue();
		}
		
		/**
		 * Get count of reads not contained in introns
		 * @return
		 */
		public Optional<Long> getCountExons() {
			try {
				return Optional.of(Long.valueOf(getGeneSpanCount() - getIntronCountAll()));
			} catch(NoSuchElementException e) {
				return Optional.empty();
			}
		}
		
		public long getGeneSpanCount() {
			return geneSpanCount.get().longValue();
		}
		
		@SuppressWarnings("unused")
		public long getGeneSpanCountMinusAllTranscribedExons() {
			return geneSpanCountMinusAllTranscribedExons.get().longValue();
		}
		
		@SuppressWarnings("unused")
		public long getGeneSpanCountMinusAllTranscribedExonsSpliced() {
			return geneSpanCountMinusAllTranscribedExonsSpliced.get().longValue();
		}
		
		public void setIntronCountUnspliced(long intronCountUnspliced) {
			this.intronCountUnspliced = Optional.of(Long.valueOf(intronCountUnspliced));
		}
		
		public void setIntronCountAll(long intronCountAll) {
			this.intronCountAll = Optional.of(Long.valueOf(intronCountAll));
		}
		
		public void setGeneSpanCountMinusAllTranscribedExons(long geneSpanCountMinusAllTranscribedExons) {
			this.geneSpanCountMinusAllTranscribedExons = Optional.of(Long.valueOf(geneSpanCountMinusAllTranscribedExons));
		}
		
		public void setGeneSpanCountMinusAllTranscribedExonsSpliced(long geneSpanCountMinusAllTranscribedExonsSpliced) {
			this.geneSpanCountMinusAllTranscribedExonsSpliced = Optional.of(Long.valueOf(geneSpanCountMinusAllTranscribedExonsSpliced));
		}
		
		public void setGeneSpanCount(long geneSpanCount) {
			this.geneSpanCount = Optional.of(Long.valueOf(geneSpanCount));
		}
		
				
		public Optional<Float> getProportionOfReads(Optional<Long> count) {
			if(count.isPresent() && geneSpanCount.isPresent()) {
				long g = geneSpanCount.get();
				long e = count.get();
				if(g != 0) return Optional.of(Float.valueOf((float)e / g));
				else return Optional.empty();
			}
			throw new NoSuchElementException();
		}
		
		public static final String columnHeaders =
				"sample_id\t"
				+ "transcript_id\t"
				+"gene_span_length\t"
				+ "num_exons\t"
				+ "total_length_of_exons\t"
				+ "total_count_over_gene_span\t"
				+ "total_count_over_exons\t"
				+ "total_count_over_introns\t"
				+ "total_count_over_introns_unspliced\t"
				+ "total_count_over_gene_span_minus_all_exons\t"
				+ "total_count_over_gene_span_minus_all_exons_spliced\t"
				+ "proportion_exons\t"
				+ "proportion_introns\t"
				+ "proportion_introns_unspliced\t";
	

		public String toString() {
			String rtrn = sampleId + "\t";
			rtrn += id + "\t";
			if(geneSpan.isPresent()) rtrn += geneSpan.get() + "\t";
			else rtrn += "-\t";
			if(numExons.isPresent()) rtrn += numExons.get() + "\t";
			else rtrn += "-\t";
			if(totalExonSize.isPresent()) rtrn += totalExonSize.get() + "\t";
			else rtrn += "-\t";
			if(geneSpanCount.isPresent()) rtrn += geneSpanCount.get() + "\t";
			else rtrn += "-\t";
			if(getCountExons().isPresent()) rtrn += getCountExons().get() + "\t";
			else rtrn += "-\t";
			if(intronCountAll.isPresent()) rtrn += intronCountAll.get() + "\t";
			else rtrn += "-\t";
			if(intronCountUnspliced.isPresent()) rtrn += intronCountUnspliced.get() + "\t";
			else rtrn += "-\t";
			if(geneSpanCountMinusAllTranscribedExons.isPresent()) rtrn += geneSpanCountMinusAllTranscribedExons.get() + "\t";
			else rtrn += "-\t";
			if(geneSpanCountMinusAllTranscribedExonsSpliced.isPresent()) rtrn += geneSpanCountMinusAllTranscribedExonsSpliced.get() + "\t";
			else rtrn += "-\t";
			try {
				rtrn += getProportionOfReads(getCountExons()).get() + "\t";
				rtrn += getProportionOfReads(intronCountAll).get() + "\t";
				rtrn += getProportionOfReads(intronCountUnspliced).get() + "\t";
			} catch(NoSuchElementException e) {
				rtrn += "-\t-\t";
			}
			return rtrn;
		}
		

		private Optional<Integer> numExons;
		private Optional<Integer> totalExonSize;
		private Optional<Integer> geneSpan;
		private Optional<Long> intronCountUnspliced;
		private Optional<Long> intronCountAll;
		private Optional<Long> geneSpanCount;
		private Optional<Long> geneSpanCountMinusAllTranscribedExons;
		private Optional<Long> geneSpanCountMinusAllTranscribedExonsSpliced;
		
	}
	
	private static Logger logger = Logger.getLogger(RatioIntronGeneSpan.class.getName());
	
	/**
	 * Get count of reads that are FULLY CONTAINED in introns
	 * @param transcript Transcript
	 * @param unsplicedOnly Unspliced reads only
	 * @return Count of reads that are FULLY CONTAINED in introns.
	 * Spliced reads that are nevertheless contained in the intron are NOT COUNTED.
	 */
	private long getTotalCountOverIntrons(BlockedAnnotation transcript, boolean unsplicedOnly) {
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
				.mapToLong(block -> getCount(chr, block.getReferenceStartPosition(), block.getReferenceEndPosition(), 
						unsplicedOnly ? Collections.singleton(isSpliced) : Collections.emptySet()))
				.sum();
	}
	
	

	private long getCount(String chr, int start, int end,
			Collection<Predicate<SAMRecord>> removeIfTrue) {
		SAMRecordIterator iter = samReader.queryContained(chr, start, end);
		long rtrn = 0;
		while(iter.hasNext()) {
			SAMRecord record = iter.next();
			boolean passesFilters = true;
			for(Predicate<SAMRecord> predicate : removeIfTrue) {
				if(predicate.test(record)) {
					passesFilters = false;
					break;
				}
			}
			if(passesFilters) rtrn++;
		}
		iter.close();
		return rtrn;
	}
		
	private static final Predicate<SAMRecord> isSpliced = (record -> {
		Iterator<CigarElement> cigarIter = record.getCigar().iterator();
		int numBlocks = 0;
		while(cigarIter.hasNext()) {
			if(cigarIter.next().getOperator().equals(CigarOperator.M)) {numBlocks++;}
		}
		return numBlocks > 1;
	});
	
	private static final Predicate<SAMRecord> isNotSpliced = (record -> !isSpliced.test(record));
	
	private long getCountOverGeneSpan(BlockedAnnotation transcript) {
		return getCount(transcript.getReferenceName(), transcript.getReferenceStartPosition(), 
				transcript.getReferenceEndPosition(), Collections.emptySet());
	}
	
	private static int getTotalSizeOfExons(BlockedAnnotation transcript) {
		return transcript.getBlockSet().stream()
				.mapToInt(block -> block.size())
				.sum();
	}
	
	private static Predicate<SAMRecord> overlapsAnyExon(AnnotationCollection<? extends Annotation> transcripts) {
		return (record -> {
			List<AlignmentBlock> blocks = record.getAlignmentBlocks();
			for(AlignmentBlock block : blocks) {
				SingleInterval blockInterval = new SingleInterval(record.getReferenceName(),
						block.getReferenceStart(), block.getReferenceStart() + block.getLength(), Strand.BOTH);
				if(transcripts.overlaps(blockInterval)) return true;
			}
			return false;
		});
	}
	
	private long getGeneSpanCountMinusExonOverlappers(Annotation transcript) {
		return getCount(transcript.getReferenceName(), transcript.getReferenceStartPosition(),
				transcript.getReferenceEndPosition(), 
				Collections.singleton(overlapsAnyExon(transcriptome)));
	}

	private long getGeneSpanCountMinusExonOverlappersSpliced(Annotation transcript) {
		Collection<Predicate<SAMRecord>> predicates = new ArrayList<Predicate<SAMRecord>>();
		predicates.add(overlapsAnyExon(transcriptome));
		predicates.add(isNotSpliced);
		return getCount(transcript.getReferenceName(), transcript.getReferenceStartPosition(),
				transcript.getReferenceEndPosition(), 
				predicates);
	}
	
	private TranscriptStats getTranscriptStats(BlockedAnnotation transcript) {
		TranscriptStats rtrn = new TranscriptStats(transcript.getName());
		rtrn.setIntronCountAll(getTotalCountOverIntrons(transcript, false));
		rtrn.setIntronCountUnspliced(getTotalCountOverIntrons(transcript, true));
		rtrn.setGeneSpan(transcript.getReferenceEndPosition() - transcript.getReferenceStartPosition());
		rtrn.setGeneSpanCount(getCountOverGeneSpan(transcript));
		rtrn.setNumExons(transcript.getNumberOfBlocks());
		rtrn.setTotalExonSize(getTotalSizeOfExons(transcript));
		rtrn.setGeneSpanCountMinusAllTranscribedExons(getGeneSpanCountMinusExonOverlappers(transcript));
		rtrn.setGeneSpanCountMinusAllTranscribedExonsSpliced(getGeneSpanCountMinusExonOverlappersSpliced(transcript));
		return rtrn;
	}
	
	private void writeTranscriptStats(BlockedAnnotation transcript, FileWriter writer) {
		try {
			TranscriptStats stats = getTranscriptStats(transcript);
			if(stats.getGeneSpanCount() == 0) return;
			writer.write(stats.toString() + "\n");
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	private void writeTranscriptStats(FileWriter writer) {
		CountLogger countLogger = new CountLogger(transcriptome.getNumAnnotations(), 10);
		CloseableIterator<? extends BlockedAnnotation> iter = transcriptome.sortedIterator();
		while(iter.hasNext()) {
			countLogger.advance();
			writeTranscriptStats(iter.next(), writer);
		}
		iter.close();
	}
	
	private void writeTranscriptStats(File outputTable) {
		try {
			FileWriter writer = new FileWriter(outputTable);
			writer.write(TranscriptStats.columnHeaders + "\n");
			writeTranscriptStats(writer);
			writer.close();
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
		p.addStringArg("-s", "Sample ID", true);
		p.parse(args);
		File bedFile = new File(p.getStringArg("-t"));
		CoordinateSpace coordSpace = new CoordinateSpace(p.getStringArg("-c"));
		File bamFile = new File(p.getStringArg("-b"));
		File output = new File(p.getStringArg("-o"));
		String sampleId = p.getStringArg("-s");
		
		RatioIntronGeneSpan r = new RatioIntronGeneSpan(bamFile, bedFile, coordSpace, sampleId);
		r.writeTranscriptStats(output);
		try {
			r.samReader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		logger.info("");
		logger.info("All done");
		
	}
	
}
