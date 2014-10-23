package motif;

import guttmanlab.core.util.CommandLineParser;
import guttmanlab.core.util.StringParser;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import broad.core.motif.SearchException;
import broad.core.motif.SequenceMotif;
import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import broad.core.sequence.SequenceRegion;

/**
 * @author prussell
 *
 */
public class SequenceMotifSearch {
	
	private Collection<SequenceMotif> motifs;
	private Collection<Sequence> sequences;
	private String mName;
	private static Logger logger = Logger.getLogger(SequenceMotifSearch.class.getName());
	
	/**
	 * @param name Experiment/sample name
	 * @param motifListFile File containing one consensus motif sequence per line
	 * @param sequenceFasta Fasta file of sequences to search
	 * @throws IOException
	 * @throws SearchException
	 */
	public SequenceMotifSearch(String name, String motifListFile, String sequenceFasta) throws IOException, SearchException {
		this(name, getFileLinesAsMotifs(motifListFile), getSequencesFromFastaFile(sequenceFasta));
	}
	
	/**
	 * @param name Experiment/sample name
	 * @param motifSet Motifs
	 * @param sequenceSet Sequences to search
	 */
	public SequenceMotifSearch(String name, Collection<SequenceMotif> motifSet, Collection<Sequence> sequenceSet) {
		motifs = motifSet;
		sequences = sequenceSet;
		mName = name;
		logger.info("There are " + motifs.size() + " motifs and " + sequences.size() + " sequences to search: " + mName + ".");
	}
	
	
	/**
	 * @return Name
	 */
	public String getName() {
		return mName;
	}
	
	/**
	 * Count all motif occurrences in sequences
	 * @return The total number of motif occurrences in all sequences
	 */
	public int countAllMotifMatches() {
		return findAllMotifMatches().size();
	}
	
	/**
	 * @return All motif matches by motif
	 */
	public Map<SequenceMotif, Collection<SequenceRegion>> findMatchesByMotif() {
		logger.info("Searching for matches by motif: " + mName + "...");
		Map<SequenceMotif, Collection<SequenceRegion>> rtrn = new HashMap<SequenceMotif, Collection<SequenceRegion>>();
		for(SequenceMotif motif : motifs) {
			Collection<SequenceRegion> matches = new ArrayList<SequenceRegion>();
			for(Sequence sequence : sequences) {
				matches.addAll(motif.match(sequence));
			}
			rtrn.put(motif, matches);
		}
		return rtrn;
	}
	
	private static Collection<SequenceMotifSearch> createFromFile(String file) throws IOException, SearchException {
		logger.info("Creating from file " + file + "...");
		FileReader r = new FileReader(file);
		BufferedReader b = new BufferedReader(r);
		StringParser s = new StringParser();
		Collection<SequenceMotifSearch> rtrn = new ArrayList<SequenceMotifSearch>();
		while(b.ready()) {
			s.parse(b.readLine());
			String name = s.asString(0);
			String motifFile = s.asString(1);
			String sequenceFasta = s.asString(2);
			rtrn.add(new SequenceMotifSearch(name, motifFile, sequenceFasta));
		}
		r.close();
		b.close();
		logger.info("Done creating from file.");
		return rtrn;
	}
	
	/**
	 * Write table of match counts
	 * @param sms SequenceMotifSearch objects
	 * @param outFile Output file
	 * @throws IOException
	 */
	public static void writeCountTable(Collection<SequenceMotifSearch> sms, String outFile) throws IOException {
		logger.info("Writing motif count table to file " + outFile + "...");
		Collection<SequenceMotif> allMotifs = new HashSet<SequenceMotif>();
		FileWriter w = new FileWriter(outFile);
		String header = "motif_id\tconsensus_sequence\tnum_possible_kmers";
		Map<SequenceMotifSearch, Map<SequenceMotif, Collection<SequenceRegion>>> matchesBySms = new HashMap<SequenceMotifSearch, Map<SequenceMotif, Collection<SequenceRegion>>>();
		for(SequenceMotifSearch s : sms) {
			allMotifs.addAll(s.getMotifs());
			header += "\t" + s.getName();
			matchesBySms.put(s, s.findMatchesByMotif());
		}
		w.write(header + "\n");
		for(SequenceMotif motif : allMotifs) {
			logger.info(motif.toString());
			String line = motif.getId() + "\t" + motif.getMotif().toString() + "\t" + motif.getNumPossibleKmers();
			for(SequenceMotifSearch s : sms) {
				if(matchesBySms.get(s).containsKey(motif)) {
					line += "\t" + matchesBySms.get(s).get(motif).size();
				} else {
					line += "\t-";
				}
			}
			w.write(line + "\n");
		}
		w.close();
		logger.info("Done writing file.");
	}
	
	/**
	 * @return The motifs
	 */
	public Collection<SequenceMotif> getMotifs() {
		return motifs;
	}
	
	/**
	 * @return Get all motif matches by sequence
	 */
	public Map<Sequence, Collection<SequenceRegion>> findMotifMatchesBySequence() {
		logger.info("Searching for matches by sequence " + mName + "...");
		Map<Sequence, Collection<SequenceRegion>> rtrn = new HashMap<Sequence, Collection<SequenceRegion>>();
		for(Sequence sequence : sequences) {
			Collection<SequenceRegion> matches = new ArrayList<SequenceRegion>();
			for(SequenceMotif motif : motifs) {
				matches.addAll(motif.match(sequence));
				for(SequenceRegion match : matches) {
					logger.debug("MATCH\t" + motif.toString() + "\t" + sequence.getId() + ":" + match.getRegionStart() + "-" + match.getRegionEnd());
				}
			}
			rtrn.put(sequence, matches);
		}
		return rtrn;
	}
	
	/**
	 * @return The number of sequences with at least one motif match
	 */
	public int numSequencesWithMatch() {
		logger.info("Counting sequences with at least one match " + mName + "...");
		Map<Sequence, Collection<SequenceRegion>> allMatches = findMotifMatchesBySequence();
		int rtrn = 0;
		for(Sequence seq : allMatches.keySet()) {
			if(allMatches.get(seq).size() > 0) {
				rtrn++;
			}
		}
		return rtrn;
	}
	
	/**
	 * @return All motif occurrences in sequences
	 */
	public Collection<SequenceRegion> findAllMotifMatches() {
		logger.info("Searching for all matches " + mName + "...");
		Collection<SequenceRegion> rtrn = new ArrayList<SequenceRegion>();
		Map<Sequence, Collection<SequenceRegion>> allMatches = findMotifMatchesBySequence();
		for(Sequence sequence : allMatches.keySet()) {
			rtrn.addAll(allMatches.get(sequence));
		}
		return rtrn;
	}
	
	private static Collection<Sequence> getSequencesFromFastaFile(String fasta) throws IOException {
		FastaSequenceIO fsio = new FastaSequenceIO(fasta);
		return fsio.loadAll();
	}
	
	private static Collection<SequenceMotif> getFileLinesAsMotifs(String file) throws IOException, SearchException {
		Collection<SequenceMotif> rtrn = new ArrayList<SequenceMotif>();
		Collection<String> asStrings = getLinesFromFile(file);
		StringParser s = new StringParser();
		for(String str : asStrings) {
			s.parse(str);
			if(s.getFieldCount() == 0) continue;
			String name = s.asString(0);
			String motif = s.asString(1);
			rtrn.add(new SequenceMotif(motif, name));
		}
		return rtrn;
	}
	
	private static Collection<String> getLinesFromFile(String file) throws IOException {
		Collection<String> rtrn = new ArrayList<String>();
		FileReader r = new FileReader(file);
		BufferedReader b = new BufferedReader(r);
		while(b.ready()) {
			rtrn.add(b.readLine());
		}
		r.close();
		b.close();
		return rtrn;
	}
	
	/**
	 * @param args
	 * @throws SearchException 
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException, SearchException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-s", "Fasta file of sequences", false, null);
		p.addStringArg("-i", "Input table of sequences. Each line: <name>   <motif_list_file>   <sequence_fasta_file>", false, null);
		p.addStringArg("-o", "Output table of match counts", false, null);
		p.addStringArg("-m", "File containing one motif per line. Each line: <identifier>   <consensus_sequence>", false, null);
		p.addBooleanArg("-d", "Debug logging", false, false);
		p.addStringArg("-n", "Single identifier", false, null);
		p.parse(args);
		boolean debug = p.getBooleanArg("-d");
		if(debug) {
			logger.setLevel(Level.DEBUG);
		}
		String sequenceFasta = p.getStringArg("-s");
		String motifFile = p.getStringArg("-m");
		String name = p.getStringArg("-n");
		String inputTable = p.getStringArg("-i");
		String outputTable = p.getStringArg("-o");

		if(sequenceFasta != null && motifFile != null && name != null) {
			SequenceMotifSearch sms = new SequenceMotifSearch(name, motifFile, sequenceFasta);
			int numMatches = sms.countAllMotifMatches();
			int numSeqsWithMatch = sms.numSequencesWithMatch();
			logger.info("There are " + numMatches + " total matches. " + numSeqsWithMatch + " sequences have at least one match.");
		}
		
		if(inputTable != null && outputTable != null) {
			writeCountTable(createFromFile(inputTable), outputTable);
		}
		
	}

}
