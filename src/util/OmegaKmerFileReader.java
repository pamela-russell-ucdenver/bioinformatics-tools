/**
 * 
 */
package util;

import guttmanlab.core.util.StringParser;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Gene;

import org.apache.log4j.Logger;

import broad.core.annotation.BasicLightweightAnnotation;
import broad.core.annotation.LightweightGenomicAnnotation;
import broad.core.math.Statistics;
import broad.core.siphy.EvolutionaryModel.OmegaFit;

/**
 * @author prussell
 *
 */
public class OmegaKmerFileReader {
	
	private int k;
	private Map<String, Map<Integer, OmegaFit>> scores;
	@SuppressWarnings("javadoc")
	public static Logger logger = Logger.getLogger(OmegaKmerFileReader.class.getName());
	private StringParser stringParser;
	private int chunkSize;
	private String dir;
	private Collection<String> omegaFilesLoaded;
	private boolean useCache;

	/**
	 * @param omegaDirectory Parent directory containing one subdirectory of omega files per chromosome
	 * @param kmerSize Kmer size
	 * @param fileChunkSize Omega file chunk size
	 */
	public OmegaKmerFileReader(String omegaDirectory, int kmerSize, int fileChunkSize) {
		this(omegaDirectory, kmerSize, fileChunkSize, true);
	}
	
	/**
	 * @param omegaDirectory Parent directory containing one subdirectory of omega files per chromosome
	 * @param kmerSize Kmer size
	 * @param fileChunkSize Omega file chunk size
	 * @param useScoreCache Cache scores when reading from file so each file is only loaded once
	 */
	public OmegaKmerFileReader(String omegaDirectory, int kmerSize, int fileChunkSize, boolean useScoreCache) {
		k = kmerSize;
		dir = omegaDirectory;
		chunkSize = fileChunkSize;
		stringParser = new StringParser();
		scores = new TreeMap<String, Map<Integer, OmegaFit>>();
		omegaFilesLoaded = new TreeSet<String>();
		useCache = useScoreCache;
		if(useCache) {
			logger.warn("Using cache for omega scores");
		} else {
			logger.warn("Not using cache for omega scores");
		}
	}
	
	private OmegaFit getFromFile(String omegaFile, String chr, int pos) throws IOException {
		if(!(pos >= getChunkBounds(omegaFile)[0] && pos <= getChunkBounds(omegaFile)[1])) {
			throw new IllegalArgumentException("Position " + pos + " is not covered by file " + omegaFile + ".");
		}
		FileReader r = new FileReader(omegaFile);
		BufferedReader b = new BufferedReader(r);
		while(b.ready()) {
			String line = b.readLine();
			stringParser.parse(line);
			if(stringParser.asInt(0) == pos) {
				r.close();
				b.close();
				return createFromDataLine(line, chr);
			}
		}
		r.close();
		b.close();
		// Position is in the interval covered by the file but does not have a score in the file
		return null;
	}
	
	private OmegaFit createFromDataLine(String dataLine, String chr) {
		stringParser.parse(dataLine);
		int pos = stringParser.asInt(0);
		LightweightGenomicAnnotation region = new BasicLightweightAnnotation(chr, pos, pos + k);
		double omega = stringParser.asDouble(1);
		double treeLength = stringParser.asDouble(2);
		double logOddsScore = stringParser.asDouble(3);
		double pval = stringParser.asDouble(4);
		OmegaFit of = new OmegaFit(region, omega, logOddsScore, pval, treeLength);			
		return of;
	}
	
	/**
	 * Get all omega scores from omega file
	 * @param omegaFile File
	 * @return Map of position to omega score object
	 * @throws IOException
	 */
	private Map<Integer, OmegaFit> getFromFile(String omegaFile) throws IOException {
		Map<Integer, OmegaFit> rtrn = new TreeMap<Integer, OmegaFit>();
		File file = new File(omegaFile);
		String chr = getChrName(omegaFile);
		if(file.exists()) {
			logger.debug("Getting kmer scores from file " + omegaFile);
			FileReader r = new FileReader(omegaFile);
			BufferedReader b = new BufferedReader(r);
			while(b.ready()) {
				String line = b.readLine();
				stringParser.parse(line);
				int pos = stringParser.asInt(0);
				rtrn.put(Integer.valueOf(pos), createFromDataLine(line, chr));
			}
			r.close();
			b.close();
		}
		logger.debug("File contained " + rtrn.size() + " scores.");
		// Set score to null for positions in the interval not reported in the file
		int numAdded = 0;
		for(int pos = getChunkBounds(omegaFile)[0]; pos <= getChunkBounds(omegaFile)[1]; pos++) {
			if(!rtrn.containsKey(Integer.valueOf(pos))) {
				rtrn.put(Integer.valueOf(pos), null);
				numAdded++;
			}
		}
		logger.debug("Added null score for " + numAdded + " positions not covered in file.");
		return rtrn;
	}
	
	/**
	 * Read an omega file and store the scores
	 * @param omegaFile Omega file
	 * @throws IOException
	 */
	private void loadFile(String omegaFile) throws IOException {
		if(omegaFilesLoaded.contains(omegaFile)) {
			throw new IllegalArgumentException("File " + omegaFile + " has already been loaded.");
		}
		String chr = getChrName(omegaFile);
		if(!scores.containsKey(chr)) {
			Map<Integer, OmegaFit> scoresThisChr = new TreeMap<Integer, OmegaFit>();
			scores.put(chr, scoresThisChr);
		}
		scores.get(chr).putAll(getFromFile(omegaFile));
		omegaFilesLoaded.add(omegaFile);
	}
	
	/**
	 * Get chromosome name from omega file name
	 * @param omegaFileName Name of omega file with or without directory path
	 * @return The chromosome name found in the file name
	 */
	private String getChrName(String omegaFileName) {
		stringParser.parse(omegaFileName,"/");
		String fileName = stringParser.asString(stringParser.getFieldCount() - 1);
		stringParser.parse(fileName,"_");
		return stringParser.asString(0);
	}
	
	/**
	 * Get the subdirectory for the chromosome
	 * @param chr Chromosome name
	 * @return Subdirectory name
	 */
	private static String getSubdir(String chr) {
		return chr.replaceAll("chr", "");
	}
	
	private int[] getChunkBounds(String omegaFile) {
		stringParser.parse(omegaFile,"/");
		String tmp1 = stringParser.asString(stringParser.getFieldCount() - 1);
		stringParser.parse(tmp1, "_");
		String tmp2 = stringParser.asString(1);
		stringParser.parse(tmp2, "\\.");
		String tmp3 = stringParser.asString(0);
		stringParser.parse(tmp3, "-");
		int[] rtrn = new int[2];
		rtrn[0] = stringParser.asInt(0);
		rtrn[1] = stringParser.asInt(1);
		return rtrn;
	}
	
	/**
	 * Get first and last position of chunk containing the position
	 * @param pos The position
	 * @return First and last position of containing chunk file
	 */
	@SuppressWarnings("unused")
	private int[] getChunkBounds(int pos) {
		int[] rtrn = new int[2];
		int remainder = pos % chunkSize;
		rtrn[0] = pos - remainder;
		rtrn[1] = rtrn[0] + chunkSize - 1;
		return rtrn;
	}
	
	/**
	 * Get name of omega file containing the position
	 * @param chr Chromosome name
	 * @param pos Position
	 * @param fullPath Whether to include directory path
	 * @return Omega file name with or without path
	 */
	private String getFileName(String chr, int pos, boolean fullPath) {
		if(pos < 0) {
			throw new IllegalArgumentException("Position must be > 0");
		}
		int remainder = pos % chunkSize;
		int chunkStart = pos - remainder;
		int chunkEnd = chunkStart + chunkSize - 1;
		String fileName = chr + "_" + chunkStart + "-" + chunkEnd + ".omega";
		String path = dir + "/" + getSubdir(chr) + "/" + fileName;
		File file = new File(path);
		if(!file.exists()) {
			logger.warn("File " + path + " does not exist.");
		}
		if(fullPath) {
			return path;
		}
		return fileName;
	}
	
	/**
	 * Get the OmegaFit object for the kmer starting at the position
	 * @param chr Chromosome
	 * @param pos Beginning position of kmer
	 * @return OmegaFit for the kmer starting at the position, or null if none is provided in omega file
	 * @throws IOException
	 */
	public OmegaFit getKmerScore(String chr, int pos) throws IOException {
		if (useCache) {
			if (scores.containsKey(chr)) {
				if (scores.get(chr).containsKey(Integer.valueOf(pos))) {
					return scores.get(chr).get(Integer.valueOf(pos));
				}
			}
			loadFile(getFileName(chr, pos, true));
			if (!scores.get(chr).containsKey(Integer.valueOf(pos))) {
				logger.debug("Scores does not contain position " + chr + " "
						+ pos);
			}
			return scores.get(chr).get(Integer.valueOf(pos));
		}
		return getFromFile(getFileName(chr, pos, true), chr, pos);
	}
	
	/**
	 * Get quantile of omega scores of all kmers in the region
	 * @param region The region
	 * @param quantile The quantile e.g. 0.5
	 * @return Omega score at the specified quantile
	 * @throws IOException
	 */
	public double getOmegaQuantile(Gene region, double quantile) throws IOException {
		double rtrn = Statistics.quantile(getOmegaScoreList(region), quantile);
		//logger.debug(quantile + " omega quantile for region " + region.getName() + " is " + rtrn + ".");
		return rtrn;
	}
	
	/**
	 * Get list of omega scores of all kmers in the region
	 * @param region The region
	 * @return List of omega scores for all kmers in the region
	 * @throws IOException
	 */
	public List<Double> getOmegaScoreList(Gene region) throws IOException {
		//logger.debug("Getting list of omega scores for region " + region.getName());
		Map<Integer, OmegaFit> omegaByPos = getAllKmerScores(region);
		List<Double> rtrn = new ArrayList<Double>();
		for(OmegaFit of : omegaByPos.values()) {
			if(of == null) {
				continue;
			}
			rtrn.add(Double.valueOf(of.getOmega()));
		}
		return rtrn;
	}
	
	/**
	 * Get omega score of each kmer in region
	 * @param region The region
	 * @return Map of kmer beginning position to score object of kmer beginning at that position
	 * @throws IOException
	 */
	public Map<Integer, OmegaFit> getAllKmerScores(Gene region) throws IOException {
		logger.debug("Getting all kmer scores for annotation " + region.getName());
		Map<Integer, OmegaFit> scoreMap = new TreeMap<Integer, OmegaFit>();
		String chr = region.getChr();
		List<Integer> positions = new ArrayList<Integer>();
		// Make list of all relevant positions in the annotation
		for (Annotation block : region.getBlocks()) {
			if (block.getSize() < k) {
				continue;
			}
			for (int pos = block.getStart(); pos <= block.getEnd() - k; pos++) {
				positions.add(Integer.valueOf(pos));
			}
		}
		if (useCache) {
			for(Integer pos : positions) {
				OmegaFit of = getKmerScore(chr, pos.intValue());
				if (of == null) {
					continue;
				}
				scoreMap.put(Integer.valueOf(pos.intValue()), of);
			}
		} else {
			int firstPos = region.getStart();
			Map<Integer, OmegaFit> tmpCache = getFromFile(getFileName(chr, firstPos, true));
			logger.debug("File contains " + tmpCache.size() + " scores.");
			for(Integer pos : positions) {
				if(!tmpCache.containsKey(pos)) {
					tmpCache.putAll(getFromFile(getFileName(chr, pos.intValue(), true)));
				}
				OmegaFit of = tmpCache.get(pos);
				if (of == null) {
					continue;
				}
				scoreMap.put(Integer.valueOf(pos.intValue()), of);
			}
		}
		logger.debug("Annotation " + region.getName() + " has " + scoreMap.size() + " scores.");
		return scoreMap;
	}

	
}
