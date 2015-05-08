package fastq;

import guttmanlab.core.sequence.FastaFileIOImpl;
import guttmanlab.core.sequence.Sequence;
import guttmanlab.core.util.CommandLineParser;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

import org.apache.log4j.Logger;

import nextgen.core.alignment.SmithWatermanAlignment;
import broad.pda.seq.fastq.FastqParser;
import broad.pda.seq.fastq.FastqSequence;

public class FastqSubsequenceMatcher {
	
	private Collection<Sequence> seqs;
	private static int MAX_MISMATCHES;
	private static Logger logger = Logger.getLogger(FastqSubsequenceMatcher.class.getName());
	private Map<String, Double> matchCounts;
	private BufferedWriter matchPctWriter;
	
	/**
	 * @param querySeqsFasta Fasta file of sequences to look for in reads
	 */
	private FastqSubsequenceMatcher(String querySeqsFasta) {
		FastaFileIOImpl fastaIO = new FastaFileIOImpl();
		seqs = fastaIO.readFromFile(querySeqsFasta);
		matchCounts = new HashMap<String, Double>();
		for(Sequence seq : seqs) {
			String seqStr = seq.getName();
			matchCounts.put(seqStr, Double.valueOf(0));
		}
	}
	
	/**
	 * Get list of names of sequences that are contained in the read
	 * @param record Fastq read
	 * @return Comma separated list of sequences contained in the read, or empty string if none
	 */
	private String getMatchList(FastqSequence record) {
		String rtrn = "";
		String readSeq = record.getSequence().toUpperCase();
		for(Sequence seq : seqs) {
			String seqStr = seq.getSequenceBases().toUpperCase();
			if(SmithWatermanAlignment.containsFullLengthUngappedMatch(readSeq, seqStr, MAX_MISMATCHES)) {
				if(!rtrn.equals("")) {
					rtrn += ",";
				}
				rtrn += seq.getName();
				matchCounts.put(seq.getName(), Double.valueOf(matchCounts.get(seq.getName()).doubleValue()+1));
			}
		}
		return rtrn;
	}
	
	/**
	 * Write table of sequence matches for each read
	 * @param fastqFile Fastq file of reads
	 * @param outTable Output table of matches per read
	 * @param outPcts Output file of percentages of reads containing sequence, updated regularly
	 * @throws IOException
	 */
	private void getMatchesAndWrite(String fastqFile, String outTable, String outPcts) throws IOException {
		logger.info("");
		logger.info("Writing table of matches to " + outTable + ".");
		logger.info("Writing percentages of reads with each sequence to " + outPcts + ".");
		BufferedWriter tableWriter = new BufferedWriter(new FileWriter(outTable));
		matchPctWriter = new BufferedWriter(new FileWriter(outPcts));
		FastqParser iter = new FastqParser();
		iter.start(new File(fastqFile));
		int numDone = 0;
		while(iter.hasNext()) {
			numDone++;
			if(numDone % 10000 == 0) {
				logger.info("Finished " + numDone + " reads.");
				writeCounts(numDone);
			}
			FastqSequence record = iter.next();
			if(record == null) {
				continue;
			}
			String matches = getMatchList(record);
			record.removeAtSymbolFromName();
			tableWriter.write(record.getName() + "\t" + matches + "\n");
		}
		tableWriter.close();
		matchPctWriter.close();
		logger.info("Done writing table and count file.");
	}
	
	/**
	 * Write current percentages to existing buffered writer
	 * @param numDone Number of reads done
	 * @throws IOException
	 */
	private void writeCounts(int numDone) throws IOException {
		matchPctWriter.write("\nReads counted: " + numDone + "\n");
		for(String seq : matchCounts.keySet()) {
			matchPctWriter.write(seq + "\t" + matchCounts.get(seq).doubleValue()/(double)numDone + "\n");
		}
	}
	
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-fq", "Fastq file of reads", true);
		p.addIntArg("-mm", "Max mismatches", true);
		p.addStringArg("-ot", "Output table of matches per read", true);
		p.addStringArg("-op", "Output file of percentage of reads with matches, updated regularly", true);
		p.addStringArg("-qs", "Fasta file of query sequences", true);
		p.parse(args);
		String fastqFile = p.getStringArg("-fq");
		String outTable = p.getStringArg("-ot");
		String outPcts = p.getStringArg("-op");
		String querySeqsFasta = p.getStringArg("-qs");
		MAX_MISMATCHES = p.getIntArg("-mm");
		
		FastqSubsequenceMatcher fsm = new FastqSubsequenceMatcher(querySeqsFasta);
		fsm.getMatchesAndWrite(fastqFile, outTable, outPcts);
		
		logger.info("");
		logger.info("All done.");
		
	}
	
}
