/**
 * 
 */
package simulation;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;
import java.util.TreeSet;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import broad.core.datastructures.Pair;
import broad.core.parser.CommandLineParser;
import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import broad.pda.annotation.BEDFileParser;

import net.sf.picard.fastq.FastqRecord;
import net.sf.picard.fastq.FastqWriter;
import net.sf.picard.fastq.FastqWriterFactory;
import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Annotation.Strand;
import nextgen.core.annotation.Gene;
import nextgen.core.feature.GeneWindow;

/**
 * @author prussell
 *
 */
public class SyntheticReadSet {
	
	private Map<String, Collection<Gene>> genes;
	private Map<String, Sequence> genome;
	protected int fragmentLength;
	protected int readLength;
	private String qualLine;
	
	protected static Logger logger = Logger.getLogger(SyntheticReadSet.class.getName());
	
	private SyntheticReadSet(String genomeFasta, String geneBed, int fragmentSize, int readSize) throws IOException {
		
		logger.info("");
		logger.info("Creating synthetic read set...");
		
		genome = FastaSequenceIO.getChrSequencesFromFasta(genomeFasta);
		genes = BEDFileParser.loadDataByChr(new File(geneBed));
		readLength = readSize;
		fragmentLength = fragmentSize;
		logger.info("");
		logger.info("Fragment size is " + fragmentLength + ".");
		logger.info("Read length is " + readLength + ".");
		qualLine = "";
		for(int i=0; i<readSize; i++) {
			qualLine += "*";
		}
	}
	
	protected Sequence getChromosome(Gene gene) {
		String chr = gene.getChr();
		return genome.get(chr);
	}
	
	@SuppressWarnings("unused")
	private void writeFastqFiles(ReadGenerator readGenerator, String outFilePrefix) {
		writeFastqFiles(readGenerator, outFilePrefix, null);
	}
	
	private void writeFastqFiles(ReadGenerator readGenerator, String outFilePrefix, String chr) {
		
		String read1file = outFilePrefix + "_1.fq";
		String read2file = outFilePrefix + "_2.fq";
		
		logger.info("");
		logger.info("Writing reads to fastq files " + read1file + " and " + read2file + "...");
		
		FastqWriterFactory f = new FastqWriterFactory();
		FastqWriter writer1 = f.newWriter(new File(read1file));
		FastqWriter writer2 = f.newWriter(new File(read2file));
		
		int numDone = 0;
		
		Collection<String> chrsToUse = new TreeSet<String>();
		if(chr == null) {
			chrsToUse.addAll(genes.keySet());
		} else {
			chrsToUse.add(chr);
		}
		
		for(String c : chrsToUse) {
			logger.info(c);
			for(Gene gene : genes.get(c)) {
				writeFastq(readGenerator.getPairedReads(gene), writer1, writer2);
				numDone++;
				if(numDone % 100 == 0) {
					logger.info("Finished " + numDone + " genes.");
				}
			}
		}
		
		writer1.close();
		writer2.close();
	}
	
	private static void writeFastq(Collection<Pair<FastqRecord>> pairedReads, FastqWriter read1writer, FastqWriter read2writer) {
		for(Pair<FastqRecord> readPair : pairedReads) {
			read1writer.write(readPair.getValue1());
			read2writer.write(readPair.getValue2());
		}
	}
	
	
	protected Pair<FastqRecord> getPairedRead(Sequence chromosome, GeneWindow fragment, String extraID) {
		logger.debug("GETTING_READS\t" + fragment.toBED());
		int size = fragment.getSize();
		Strand orientation = fragment.getOrientation();
		String sourceAnnotName = "";
		for(Annotation a : fragment.getSourceAnnotations()) {
			sourceAnnotName += a.getName() + ":";
		}
		if(orientation.equals(Strand.UNKNOWN)) {
			throw new IllegalArgumentException("Strand must be known " + fragment.toBED());
		}
		Annotation read1 = orientation.equals(Strand.POSITIVE) ? fragment.trimGene(size - readLength, size) : fragment.trimGene(0, readLength);
		logger.debug("READ1\t" + read1.toBED());
		Annotation read2 = orientation.equals(Strand.POSITIVE) ? fragment.trimGene(0, readLength) : fragment.trimGene(size - readLength, size);
		logger.debug("READ2\t" + read2.toBED());
		Sequence read1seq = chromosome.getSubsequence(read1);
		read1seq.reverse();
		Sequence read2seq = chromosome.getSubsequence(read2);
		String seqHeaderPrefix1 = sourceAnnotName + fragment.getChr() + ":" + fragment.getStart() + "-" + fragment.getEnd() + ":" + orientation.toString();
		String seqHeaderPrefix2 = seqHeaderPrefix1;
		if(extraID != null) {
			seqHeaderPrefix1 += ":" + extraID;
			seqHeaderPrefix2 += ":" + extraID;
		}
		String qualHeaderPrefix = "";
		FastqRecord read1fastq = new FastqRecord(seqHeaderPrefix1, read1seq.getSequenceBases().toUpperCase(), qualHeaderPrefix, qualLine);
		FastqRecord read2fastq = new FastqRecord(seqHeaderPrefix2, read2seq.getSequenceBases().toUpperCase(), qualHeaderPrefix, qualLine);
		Pair<FastqRecord> rtrn = new Pair<FastqRecord>(read1fastq, read2fastq);
		return rtrn;
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		String description = "\n";
		description += "Simulate a read set from annotated transcripts.\n";
		description += "Task 1: Create one fragment beginning at each transcribed position.\n"; 
		description += "        No differences to transcript sequence.\n";
		description += "Task 2: Introduce deletions.\n";
		description += "        For each transcribed position, create fragments with a deletion\n";
		description += "        at each position of the fragment.\n";
		p.setProgramDescription(description);
		p.addStringArg("-gf", "Genome fasta file", true);
		p.addStringArg("-gb", "Gene annotation in bed format", true);
		p.addIntArg("-f", "Size of fragments to create", true);
		p.addIntArg("-r", "Length of paired reads to create from fragments", true);
		p.addStringArg("-o", "Output fastq file prefix", true);
		p.addBooleanArg("-ep", "Perform task 1 (one fragment for every position)", false, false);
		p.addBooleanArg("-epd", "Perform task 2 (one fragment for every deletion at every position)", false, false);
		p.addStringArg("-c", "Single chromosome to use", false, null);
		p.addIntArg("-d", "Deletion length for task 2", false, 1);
		p.addBooleanArg("--debug", "Debug logging", false, false);
		p.parse(args);
		String genomeFasta = p.getStringArg("-gf");
		String geneBed = p.getStringArg("-gb");
		String chr = p.getStringArg("-c");
		int fragmentSize = p.getIntArg("-f");
		int readLength = p.getIntArg("-r");
		String outFilePrefix = p.getStringArg("-o");
		boolean everyPosition = p.getBooleanArg("-ep");
		boolean everyPositionEveryDeletion = p.getBooleanArg("-epd");
		int deletionSize = p.getIntArg("-d");
		if(p.getBooleanArg("--debug")) {
			logger.setLevel(Level.DEBUG);
		}
		
		if((!everyPosition && !everyPositionEveryDeletion) || everyPosition && everyPositionEveryDeletion) {
			throw new IllegalArgumentException("Must provide exactly one of -ep, -epd.");
		}
		
		SyntheticReadSet s = new SyntheticReadSet(genomeFasta, geneBed, fragmentSize, readLength);
		if(everyPosition) s.writeFastqFiles(s.new ReadGeneratorEveryPosition(), outFilePrefix, chr);
		if(everyPositionEveryDeletion) s.writeFastqFiles(s.new ReadGeneratorEveryPositionEveryDeletion(deletionSize), outFilePrefix, chr);
		
		logger.info("");
		logger.info("All done.");
		
	}
	
	private class ReadGeneratorEveryPositionEveryDeletion implements ReadGenerator {

		private int deletionLength;
		
		public ReadGeneratorEveryPositionEveryDeletion(int deletionSize) {
			deletionLength = deletionSize;
		}
		
		@Override
		public Collection<Pair<FastqRecord>> getPairedReads(Gene gene) {
			if(deletionLength < 1 || deletionLength >= readLength) {
				throw new IllegalArgumentException("Deletion length must be between 1 and " + Integer.valueOf(readLength-1).toString() + ".");
			}
			Collection<Pair<FastqRecord>> rtrn = new ArrayList<Pair<FastqRecord>>();
			if(gene.size() < fragmentLength) {
				return rtrn;
			}
			Collection<GeneWindow> fragments = gene.getWindows(fragmentLength, 1, 0);
			Sequence chromosome = getChromosome(gene);
			for(GeneWindow fragment : fragments) {
				int fragmentSize = fragment.size();
				logger.debug("FRAGMENT\t" + fragment.toBED());
				for(int deletionStart = 0; deletionStart <= fragmentSize - deletionLength; deletionStart++) {
					Annotation deletion = fragment.copy();
					deletion.trim(deletionStart, fragmentSize - deletionStart - deletionLength);
					logger.debug("DELETION\tstart=" + deletionStart + "\t" + deletion.toBED()); 
					Annotation fragmentWithDeletion = fragment.minus(deletion); // correct
					logger.debug("FRAGMENT_WITH_DELETION\t" + fragmentWithDeletion.toBED());
					String extraID = "deletion_" + deletionStart + "_" + deletionLength;
					rtrn.add(getPairedRead(chromosome, new GeneWindow(fragmentWithDeletion), extraID));
				}
			}
			return rtrn;
		}
		
	}
	
	private class ReadGeneratorEveryPosition implements ReadGenerator {
		
		public ReadGeneratorEveryPosition() {}
		
		@Override
		public Collection<Pair<FastqRecord>> getPairedReads(Gene gene) {
			Collection<Pair<FastqRecord>> rtrn = new ArrayList<Pair<FastqRecord>>();
			if(gene.size() < fragmentLength) {
				return rtrn;
			}
			Sequence chromosome = getChromosome(gene);
			Collection<GeneWindow> fragments = gene.getWindows(fragmentLength, 1, 0);
			for(GeneWindow fragment : fragments) {
				rtrn.add(getPairedRead(chromosome, fragment, null));
			}
			return rtrn;
		}
		
	}
	
	private interface ReadGenerator {
		
		public Collection<Pair<FastqRecord>> getPairedReads(Gene gene);
		
	}
	

}
