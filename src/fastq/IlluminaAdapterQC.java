/**
 * 
 */
package fastq;

import general.CommandLineParser;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;



/**
 * @author prussell
 *
 *
 */
public class IlluminaAdapterQC {

	/*
	 * The length of the reads
	 * This is the maximum length of adapter fragments to search for
	 * Can be longer than the actual reads
	 */
	private int readLength;
	
	/**
	 * Fasta file containing sequences of all adapters and their reverse complements
	 */
	private String adapterSequenceFasta;
	
	/**
	 * The adapter sequences in a container
	 */
	private List<Sequence> adapterSequences;
	
	/**
	 * The fastq or fasta file containing first read mates
	 */
	private String read1file;
	
	/**
	 * The fastq or fasta file containing second read mates
	 */
	private String read2file;
	
	/**
	 * Output file to write QC metrics to
	 */
	private String outCountsFile;
	
	/**
	 * Output fasta file to write adapter-trimmed read1 to
	 */
	private String outTrimmedRead1File;
	
	/**
	 * Output fasta file to write adapter-trimmed read2 to
	 */	
	private String outTrimmedRead2File;
	
	/**
	 * Whether to write trimmed read1 to file
	 */
	private boolean writeTrimmedRead1File;
	
	/**
	 * Whether to write trimmed read2 to file
	 */
	private boolean writeTrimmedRead2File;
	
	/**
	 * The number of reads
	 */
	private int numReads;
	
	/**
	 * All subsequences of the adapters up to the read length
	 */
	private Map<Integer, TreeSet<String>> adapterSubsequences;
	
	/**
	 * Whether the input reads are in fastq format (or fasta)
	 */
	private boolean inputIsFastq;
	
	/**
	 * Counts of amount of adapter in reads
	 * Row is adapter content of read1, column is adapter content of read2, value is number of read pairs with those amounts of adapter
	 */
	private int[][] counts;
	
	/**
	 * File stream to write trimmed read1 to
	 */
	private FileWriter trimmedRead1Writer;
	
	/**
	 * File stream to write trimmed read2 to
	 */
	private FileWriter trimmedRead2Writer;
	
	/**
	 * Constructor without writing trimmed reads to fasta file
	 * @param readlength a cap on the read length
	 * @param adapterSequenceFile the path of a fasta file containing adapter sequences and their reverse complements
	 * @param read1File the path of a fastq or fasta file containing the first read mates
	 * @param read2File the path of a fastq or fasta file containing the second read mates
	 * @param fastq whether the input reads are fastq (or fasta)
	 * @param outfile the output file to write QC metrics to
	 * @throws IOException
	 */
	public IlluminaAdapterQC(int readlength, String adapterSequenceFile, String read1File, String read2File, boolean fastq, String outfile) throws IOException {
		this(readlength, adapterSequenceFile, read1File, read2File, fastq, outfile, null, null);
	}
	
	/**
	 * Constructor with writing trimmed reads to fasta file
	 * @param readlength a cap on the read length
	 * @param adapterSequenceFile the path of a fasta file containing adapter sequences and their reverse complements
	 * @param read1File the path of a fastq or fasta file containing the first read mates
	 * @param read2File the path of a fastq or fasta file containing the second read mates
	 * @param fastq whether the input reads are fastq (or fasta)
	 * @param outfile the output file to write QC metrics to
	 * @param outRead1file the output fasta file to write adapter-trimmed read1 to
	 * @param outRead2file the output fasta file to write adapter-trimmed read2 to
	 * @throws IOException
	 */
	public IlluminaAdapterQC(int readlength, String adapterSequenceFile, String read1File, String read2File, boolean fastq, String outfile, String outRead1file, String outRead2file) throws IOException {
		
		// set attributes
		this.readLength = readlength;
		this.counts = new int[readlength+1][readlength+1];
		this.outCountsFile = outfile;
		this.outTrimmedRead1File = outRead1file;
		this.outTrimmedRead2File = outRead2file;
		this.numReads = 0;
		this.inputIsFastq = fastq;
		
		// decide whether to write trimmed reads
		this.writeTrimmedRead1File = this.outTrimmedRead1File != null;
		this.writeTrimmedRead2File = this.outTrimmedRead2File != null;
		
		if(this.writeTrimmedRead1File) this.trimmedRead1Writer = new FileWriter(this.outTrimmedRead1File);
		if(this.writeTrimmedRead2File) this.trimmedRead2Writer = new FileWriter(this.outTrimmedRead2File);
		
		// initialize adapter counts
		for(int i=0; i<this.readLength+1; i++) {
			for(int j=0; j<this.readLength+1; j++) {
				this.counts[i][j] = 0;
			}
		}
		
		// read adapter sequences
		this.adapterSequenceFasta = adapterSequenceFile;
		FastaSequenceIO fsio = new FastaSequenceIO(this.adapterSequenceFasta);
		this.adapterSequences = fsio.loadAll();
		
		// initialize adapter subsequences
		this.adapterSubsequences = new TreeMap<Integer, TreeSet<String>>();
		
		// store adapter subsequences
		for(int i=1; i<this.readLength+1; i++) {
			TreeSet<String> subseqs = new TreeSet<String>();
			for(Sequence seq : this.adapterSequences) {
				subseqs.add(seq.getSequenceBases().substring(0,i));
			}
			this.adapterSubsequences.put(Integer.valueOf(i),subseqs);
		}
		
		this.read1file = read1File;
		this.read2file = read2File;
	}
	
	/**
	 * Whether to write trimmed read1 to fasta file
	 * @return Whether to write trimmed read1 to fasta file
	 */
	public boolean writeRead1() {return this.writeTrimmedRead1File;}
	
	/**
	 * Whether to write trimmed read2 to fasta file
	 * @return Whether to write trimmed read2 to fasta file
	 */	
	public boolean writeRead2() {return this.writeTrimmedRead2File;}
	
	/**
	 * Get the amount of adapter present at the end of a read
	 * @param read the read sequence
	 * @return the length of longest adapter prefix present as a suffix of the read sequence
	 */
	private int getAdapterSubsequenceLength(String read) {

		boolean found = false;
		
		for(int i = this.readLength; i > 0; i--) {
			if(read.length() < i) continue;
			if(found) break;
			for(String subseq : this.adapterSubsequences.get(Integer.valueOf(i))) {
				if(read.substring(read.length()-i,read.length()).equals(subseq)) {
					return i;
				}
			}
		}

		return 0;
		
	}

	
	private void analyzeAdapters() throws IOException {

		// reset counts
		for(int i=0; i<this.readLength+1; i++) {
			for(int j=0; j<this.readLength+1; j++) {
				this.counts[i][j] = 0;
			}
		}
		this.numReads = 0;
	
		FileReader f1 = new FileReader(this.read1file);
		FileReader f2 = new FileReader(this.read2file);
		BufferedReader b1 = new BufferedReader(f1);
		BufferedReader b2 = new BufferedReader(f2);
		
		int lineNumber = 0;
		
		while(b1.ready() && b2.ready()) {
			
			lineNumber++;
			String read1 = b1.readLine();
			String read2 = b2.readLine();
			
			// skip lines not containing read sequence
			if(this.inputIsFastq && lineNumber % 4 != 2) continue;
			if(!this.inputIsFastq && lineNumber % 2 != 0) continue;
			
			this.numReads++;
			
			// get amount of adapter in each read sequence
			int read1AdapterBases = this.getAdapterSubsequenceLength(read1);
			int read2AdapterBases = this.getAdapterSubsequenceLength(read2);
			
			// write trimmed read1 to file
			if(this.writeTrimmedRead1File) {
				this.trimmedRead1Writer.write(">read_" + this.numReads + ":1\n");
				this.trimmedRead1Writer.write(read1.substring( 0, read1.length() - read1AdapterBases ) + "\n");
			}

			// write trimmed read2 to file
			if(this.writeTrimmedRead2File) {
				this.trimmedRead2Writer.write(">read_" + this.numReads + ":2\n");
				this.trimmedRead2Writer.write(read2.substring( 0, read2.length() - read2AdapterBases ) + "\n");
			}

			// increment the count of reads with these adapter subsequence lengths
			this.counts[read1AdapterBases][read2AdapterBases]++;
			
		}
		
		if(this.writeTrimmedRead1File) this.trimmedRead1Writer.close();
		if(this.writeTrimmedRead2File) this.trimmedRead2Writer.close();
		
		b1.close();
		b2.close();
		
	}
	
	/**
	 * Print QC metrics to file
	 * Writes matrix of number of read pairs with given lengths of adapter subsequence in each read
	 * Also writes vectors of totals for read1 and read2
	 * @throws IOException
	 */
	private void printCounts() throws IOException {
		
		FileWriter w = new FileWriter(this.outCountsFile);
		
		w.write("Read1_adapter_amount\n");
		for(int i=0; i <= this.readLength; i++) {
			int sum = 0;
			for(int j=0; j <= this.readLength; j++) sum += this.counts[i][j];
			w.write(i + "\t" + sum + "\t" + (float)sum/(float)this.numReads + "\n");
		}
		
		w.write("\n\n\n");
		
		w.write("Read2_adapter_amount\n");
		for(int i=0; i <= this.readLength; i++) {
			int sum = 0;
			for(int j=0; j <= this.readLength; j++) sum += this.counts[j][i];
			w.write(i + "\t" + sum + "\t" + (float)sum/(float)this.numReads + "\n");
		}
		
		w.write("\n\n\n");
		
		w.write("\tRead2_adapter_amount->\n");
		w.write("Read1_adapter_amount\t");
		for(int i=0; i <= this.readLength; i++) {
			w.write(i + "\t");
		}
		w.write("\n");
		
		for(int i=0; i <= this.readLength; i++) {
			w.write(i + "\t");
			for(int j=0; j <= this.readLength; j++) {
				w.write(this.counts[i][j] + "\t");
			}
			w.write("\n");
		}
		
		w.write("\n\n\n");
		w.write("Read2_adapter_amount->\n");
		w.write("Read1_adapter_amount\t");
		for(int i=0; i <= this.readLength; i++) {
			w.write(i + "\t");
		}
		w.write("\n");
		
		for(int i=0; i <= this.readLength; i++) {
			w.write(i + "\t");
			for(int j=0; j <= this.readLength; j++) {
				w.write((float)this.counts[i][j]/(float)this.numReads + "\t");
			}
			w.write("\n");
		}
		
		
		w.close();
		
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {

		// Get arguments from command line
		CommandLineParser c = new CommandLineParser();
		c.addIntArg("-l", "Read length", true);
		c.addStringArg("-a", "Fasta file of adapter sequences", true);
		c.addStringArg("-1", "Fastq or fasta file of read1", true);
		c.addStringArg("-2", "Fastq or fasta file of read2", true);
		c.addStringArg("-o", "Output file for adapter counts", true);
		c.addStringArg("-of1", "Fasta file to write trimmed read1 to", false, null);
		c.addStringArg("-of2", "Fasta file to write trimmed read2 to", false, null);
		c.addBooleanArg("-fa", "Input reads are in fasta format", false, false);
		c.parse(args);

		int readLength = c.getIntArg("-l");
		String adapters = c.getStringArg("-a");
		String read1 = c.getStringArg("-1");
		String read2 = c.getStringArg("-2");
		String outfile = c.getStringArg("-o");
		String outreads1 = c.getStringArg("-of1");
		String outreads2 = c.getStringArg("-of2");
		boolean fasta = c.getBooleanArg("-fa");
		
		IlluminaAdapterQC qc = new IlluminaAdapterQC(readLength,adapters,read1,read2,!fasta,outfile,outreads1,outreads2);
		
		qc.analyzeAdapters();
		qc.printCounts();
		
	}

}
