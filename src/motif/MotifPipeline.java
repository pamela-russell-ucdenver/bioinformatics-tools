/**
 * 
 */
package motif;

import general.CommandLineParser;
import general.StringParser;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;

import org.apache.log4j.Logger;
import org.ggf.drmaa.DrmaaException;
import org.ggf.drmaa.Session;

import pipeline.Job;
import pipeline.JobUtils;
import pipeline.LSFJob;
import pipeline.OGSJob;
import pipeline.Scheduler;

import annotation.WindowWriter;
import bed.BedFileCollapseOverlappers;
import bed.BedFileFilter;
import broad.core.math.EmpiricalDistribution;
import broad.pda.annotation.BEDFileParser;

import nextgen.core.annotation.Gene;
import nextgen.core.pipeline.util.OGSUtils;
import nextgen.core.programs.FastaAnnotationExtractor;

/**
 * @author prussell
 *
 */
public class MotifPipeline {
	
	private String genomeFasta;
	private int medianFeatureSize;
	
	private String featureBed;
	protected Map<String, Collection<Gene>> features;
	private Map<String, Collection<Gene>> genes;
	private String collapsedFeatureBed;
	private String collapsedFeatureFasta;
	private String geneBed;
	
	private Session session;
	
	private String genesMinusCollapsedFeaturesBed;
	private String genesMinusCollapsedFeaturesWindowBed;
	private String genesMinusCollapsedFeaturesWindowFasta;

	private static String IND_GENE_INPUT_DIRECTORY = "individual_gene_input_files";
	private static String BSUB_DIR = "bsub_output";
	private static String DREME_DIR = "dreme_output";
	private static String FIMO_DIR = "fimo_output";
	private static String ALL_FEATURES_JOB_DESCRIPTION = "all_features";
	private static double ESTIMATE_FIMO_ALPHA_INDIVIDUAL_GENE = 0.25;
	private static double DEFAULT_FIMO_QVAL_THRESHOLD = 0.05;

	private static Logger logger = Logger.getLogger(MotifPipeline.class.getName());
	private ArrayList<String> dremeJobDescriptions;
	private ArrayList<String> fimoJobDescriptions;

	/**
	 * @param geneBedFile Bed file of genes
	 * @param featureBedFile Bed file of features
	 * @param genomeFastaFile Genome fasta file
	 * @throws IOException
	 */
	public MotifPipeline(String geneBedFile, String featureBedFile, String genomeFastaFile, Session drmaaSession) throws IOException {

		session = drmaaSession;
		
		File d = new File(BSUB_DIR);
		@SuppressWarnings("unused")
		boolean madeDir = d.mkdir();
		File d2 = new File(IND_GENE_INPUT_DIRECTORY);
		@SuppressWarnings("unused")
		boolean madeDir2 = d2.mkdir();
		File d3 = new File(DREME_DIR);
		@SuppressWarnings("unused")
		boolean madeDir3 = d3.mkdir();
		File d4 = new File(FIMO_DIR);
		@SuppressWarnings("unused")
		boolean madeDir4 = d4.mkdir();
		dremeJobDescriptions = new ArrayList<String>();
		fimoJobDescriptions = new ArrayList<String>();

		// Load genes, features and genome
		logger.info("");
		logger.info("Loading features and genome...");
		genes = BEDFileParser.loadDataByChr(new File(geneBedFile));
		features = BEDFileParser.loadDataByChr(new File(featureBedFile));
		genomeFasta = genomeFastaFile;
		
		// Shorten feature names if necessary
		StringParser s = new StringParser();
		for(String chr : features.keySet()) {
			for(Gene feature : features.get(chr)) {
				String longName = feature.getName();
				if(longName.length() < 150) continue;
				s.parse(longName, ":");
				String newName = null;
				if(s.getFieldCount() > 1) {
					String longNameFirstPart = s.asString(0);
					String longNameShortenFirstPart = longNameFirstPart.substring(0,100);
					newName = longName.replaceAll(longNameFirstPart, longNameShortenFirstPart);
				} else {
					newName = longName.substring(0, 100);
				}
				logger.warn("Renaming feature " + longName + " to " + newName);
				feature.setName(newName);			
			}
		}
		
		// Write new bed file of features with shortened names
		String shortenedNameFeatureBed = featureBedFile.replaceAll(".bed", "") + "_shortened_names.bed";
		FileWriter fw = new FileWriter(shortenedNameFeatureBed);
		for(String chr : features.keySet()) {
			for(Gene feature : features.get(chr)) {
				fw.write(feature.toBED() + "\n");
			}
		}
		fw.close();
		
		// Replace illegal characters
		featureBed = shortenedNameFeatureBed.replaceAll(".bed", "") + "_validChars.bed";
		replaceCharacter(shortenedNameFeatureBed, featureBed, ':', '_');		
		
		// Shorten gene names if necessary
		for(String chr : genes.keySet()) {
			for(Gene gene : genes.get(chr)) {
				String longName = gene.getName();
				if(longName.length() >= 170) {
					String newName = longName.substring(0, 170);
					logger.warn("Renaming gene " + longName + " to " + newName);
					gene.setName(newName);			
				}
			}
		}
		
		// Write new bed file of genes with shortened names
		geneBed = geneBedFile.replaceAll(".bed", "") + "_shortNames.bed";
		FileWriter w = new FileWriter(geneBed);
		for(String chr : genes.keySet()) {
			for(Gene gene : genes.get(chr)) {
				w.write(gene.toBED() + "\n");
			}
		}
		w.close();
		
		// Get median feature size
		logger.info("");
		logger.info("Getting median feature size...");
		ArrayList<Double> featureSizes = new ArrayList<Double>();
		for(String chr : features.keySet()) {
			for(Gene feature : features.get(chr)) {
				featureSizes.add(Double.valueOf(feature.getSize()));
			}
		}
		EmpiricalDistribution featureSizeDist = new EmpiricalDistribution(featureSizes);
		medianFeatureSize = (int) featureSizeDist.getMedianOfAllDataValues();
		logger.info("Median feature size is " + medianFeatureSize + ".");
		
		// Collapse overlapping features
		logger.info("");
		logger.info("Collapsing overlapping features...");
		collapsedFeatureBed = featureBed.replaceAll(".bed", "") + "_collapsed.bed";
		File f1 = new File(collapsedFeatureBed);
		if(f1.exists()) {
			logger.warn("File " + f1 + " already exists. Not remaking file.");
		} else {
			BedFileCollapseOverlappers.collapseOverlappersAndWrite(featureBed, collapsedFeatureBed, false);
		}
		
		// Write fasta file of collapsed features
		logger.info("");
		logger.info("Writing fasta file of collapsed features...");
		collapsedFeatureFasta = collapsedFeatureBed.replaceAll(".bed", "") + ".fa";
		File f2 = new File(collapsedFeatureFasta);
		if(f2.exists()) {
			logger.warn("File " + f2 + " already exists. Not remaking file.");
		} else {
			FastaAnnotationExtractor fa = new FastaAnnotationExtractor(genomeFastaFile, collapsedFeatureBed);
			fa.writeFasta(collapsedFeatureFasta);
		}
		
		// Make bed file of genes minus collapsed features
		logger.info("");
		logger.info("Writing bed file of genes minus collapsed features...");
		genesMinusCollapsedFeaturesBed = geneBedFile.replaceAll(".bed", "") + "_minusFeatures.bed";
		File f3 = new File(genesMinusCollapsedFeaturesBed);
		if(f3.exists()) {
			logger.warn("File " + f3 + " already exists. Not remaking file.");
		} else {
			BedFileFilter bff = new BedFileFilter(geneBedFile);
			bff.subtractSegmentsFromGenes(collapsedFeatureBed, genesMinusCollapsedFeaturesBed, false);
		}
		
		// Write windows of genes minus peaks
		logger.info("");
		logger.info("Writing windows of genes minus collapsed features...");
		String genesMinusCollapsedFeaturesWindowBedInvalidChars = genesMinusCollapsedFeaturesBed.replaceAll(".bed", "") + "_window_" + medianFeatureSize + ".bed";
		genesMinusCollapsedFeaturesWindowBed = genesMinusCollapsedFeaturesWindowBedInvalidChars.replaceAll(".bed", "") + "_validChars.bed";
		File f6 = new File(genesMinusCollapsedFeaturesWindowBed);
		if(f6.exists()) {
			logger.warn("File " + f6 + " already exists. Not remaking file.");
		} else {
			File f4 = new File(genesMinusCollapsedFeaturesWindowBedInvalidChars);
			if(f4.exists()) {
				logger.warn("File " + f4 + " already exists. Not remaking file.");
			} else {
				WindowWriter ww = new WindowWriter(genesMinusCollapsedFeaturesBed);
				ww.writeWindowsToFile(medianFeatureSize, medianFeatureSize, genesMinusCollapsedFeaturesWindowBedInvalidChars);
			}
			replaceCharacter(genesMinusCollapsedFeaturesWindowBedInvalidChars, genesMinusCollapsedFeaturesWindowBed, ':', '_');
		}
		
		// Write sequences of genes minus peaks (windows)
		logger.info("");
		logger.info("Writing sequences of windows of genes minus collapsed features...");
		genesMinusCollapsedFeaturesWindowFasta = genesMinusCollapsedFeaturesWindowBed.replaceAll(".bed", "") + ".fa";
		File f5 = new File(genesMinusCollapsedFeaturesWindowFasta);
		if(f5.exists()) {
			logger.warn("File " + f5 + " already exists. Not remaking file.");
		} else {
			FastaAnnotationExtractor fa2 = new FastaAnnotationExtractor(genomeFasta, genesMinusCollapsedFeaturesWindowBed);
			fa2.writeFasta(genesMinusCollapsedFeaturesWindowFasta);
		}
		
		logger.info("");
		logger.info("Done preparing sequences.");
		
	}
	
	private static void replaceCharacter(String inputFile, String outputFile, char oldChar, char newChar) throws IOException {
		FileReader r = new FileReader(inputFile);
		BufferedReader b = new BufferedReader(r);
		FileWriter w = new FileWriter(outputFile);
		while(b.ready()) {
			String line = b.readLine();
			String newline = line.replace(oldChar, newChar);
			w.write(newline + "\n");
		}
		r.close();
		b.close();
		w.close();
	}
	
	private static String getIndividualGeneBedFileName(Gene gene) {
		return IND_GENE_INPUT_DIRECTORY + "/" + gene.getName() + ".bed";
	}
	
	private static String getIndividualGeneOverlappingFeaturesBedFileName(Gene gene) {
		return IND_GENE_INPUT_DIRECTORY + "/" + gene.getName() + "_overlapping_features.bed";
	}
	
	private static String getIndividualGeneOverlappingCollapsedFeaturesBedFileName(Gene gene) {
		return IND_GENE_INPUT_DIRECTORY + "/" + gene.getName() + "_overlapping_features_shortened_names_validChars_collapsed.bed";
	}
	
	/**
	 * Make bed file of one gene and features that overlap the gene
	 * @param gene The gene
	 * @throws IOException
	 */
	private void makeIndividualGeneBedFiles(Gene gene) throws IOException {
		String bed = getIndividualGeneBedFileName(gene);
		String featureBed1 = getIndividualGeneOverlappingFeaturesBedFileName(gene);
		FileWriter w = new FileWriter(bed);
		w.write(gene.toBED() + "\n");
		w.close();
		FileWriter w2 = new FileWriter(featureBed1);
		if(!features.containsKey(gene.getChr())) {
			w2.close();
			return;
		}
		for(Gene feature : features.get(gene.getChr())) {
			if(feature.overlaps(gene)) {
				w2.write(feature.toBED() + "\n");
			}
		}
		w2.close();
	}

	/**
	 * Run DREME on features within one gene
	 * @param gene The gene
	 * @param dremeExecutable DREME executable
	 * @param drmaaSession DRMAA session or null if not using OGS
	 * @return LSF job ID
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws DrmaaException 
	 */
	private Job submitDremeOnGene(Gene gene, String dremeExecutable, String batchedMotifPipelineJar, Scheduler scheduler, Session drmaaSession) throws IOException, InterruptedException, DrmaaException {
		String cmmd = "java -jar -Xmx14g -Xms10g -Xmn8g ";
		cmmd += batchedMotifPipelineJar + " ";
		cmmd += "-gb " + getIndividualGeneBedFileName(gene) + " ";
		cmmd += "-fb " + getIndividualGeneOverlappingFeaturesBedFileName(gene) + " ";
		cmmd += "-g " + genomeFasta + " ";
		cmmd += "-d " + dremeExecutable + " ";
		cmmd += "-a " + ESTIMATE_FIMO_ALPHA_INDIVIDUAL_GENE + " ";
		cmmd += "-q " + DEFAULT_FIMO_QVAL_THRESHOLD + " ";
		String description = getJobDescription(gene);
		cmmd += "-de " + description;
		dremeJobDescriptions.add(description);
		String jobID = "dreme_" + description + "_" + System.currentTimeMillis();
		switch(scheduler) {
		case LSF:
			String bsubOut = BSUB_DIR + "/" + jobID + ".bsub";
			LSFJob job = new LSFJob(Runtime.getRuntime(), jobID, cmmd, bsubOut, "hour", 16);
			job.submit();
			return job;
		case OGS:
			OGSJob ogsjob = new OGSJob(drmaaSession, cmmd);
			ogsjob.submit();
			return ogsjob;
		default:
			throw new IllegalArgumentException("Scheduler " + scheduler.toString() + " not supported.");
		}
	}
	
	/**
	 * Get job description for gene
	 * @param gene Gene
	 * @return Description
	 */
	private static String getJobDescription(Gene gene) {
		return gene.getName();
	}
	
	/**
	 * Run Fimo on motifs found within one gene
	 * @param gene The gene
	 * @param fimoExecutable Fimo executable
	 * @param fimoOptionAlpha Expected proportion of sequences with the motif
	 * @param qvalThreshold Q value threshold
	 * @param DRMAA session or null if not using OGS
	 * @return LSF job ID
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws DrmaaException 
	 */
	private Job submitFimoOnGene(Gene gene, String fimoExecutable, double fimoOptionAlpha, double qvalThreshold, String batchedMotifPipelineJar, Scheduler scheduler, Session drmaaSession) throws IOException, InterruptedException, DrmaaException {
		String cmmd = "java -jar -Xmx3g -Xms2g -Xmn1g ";
		cmmd += batchedMotifPipelineJar + " ";
		cmmd += "-f " + fimoExecutable + " ";
		cmmd += "-gb " + getIndividualGeneBedFileName(gene) + " ";
		cmmd += "-fb " + getIndividualGeneOverlappingFeaturesBedFileName(gene) + " ";
		cmmd += "-g " + genomeFasta + " ";
		cmmd += "-a " + fimoOptionAlpha + " ";
		cmmd += "-q " + qvalThreshold + " ";
		String description = getJobDescription(gene);
		cmmd += "-de " + description;
		fimoJobDescriptions.add(description);
		String jobID = "fimo_" + description + "_" + System.currentTimeMillis();
		switch(scheduler) {
		case LSF:
			String bsubOut = BSUB_DIR + "/" + jobID + ".bsub";
			LSFJob job = new LSFJob(Runtime.getRuntime(), jobID, cmmd, bsubOut, "hour", 4);
			job.submit();
			return job;
		case OGS:
			OGSJob ogsjob = new OGSJob(drmaaSession, cmmd);
			ogsjob.submit();
			return ogsjob;
		default:
			throw new IllegalArgumentException("Scheduler " + scheduler + " not supported.");
		}
	}
	
	/**
	 * Whether the DREME run found any motifs
	 * @param jobDescription DREME job description
	 * @return True iff there the motif file contains motifs
	 */
	private static boolean dremeFoundMotif(String jobDescription) throws IOException {
		String dremeFile = getDremeDirectory(jobDescription) + "/dreme.txt";
		File f = new File(dremeFile);
		if(!f.exists()) return false;
		FileReader r = new FileReader(dremeFile);
		BufferedReader b = new BufferedReader(r);
		boolean rtrn = false;
		while(b.ready()) {
			String line = b.readLine();
			if(line.contains("MOTIF")) {
				rtrn = true;
			}
		}
		r.close();
		b.close();
		return rtrn;
	}
	
	/**
	 * Whether the FIMO job found any motif instances
	 * @param jobDescription FIMO job description
	 * @return True iff the fimo.txt file has data
	 * @throws IOException
	 */
	private static boolean fimoFoundMotif(String jobDescription) throws IOException {
		String fimoDir = getFimoDirectory(jobDescription);
		File dirFile = new File(fimoDir);
		if(!dirFile.exists()) return false;
		String fimoTxt = fimoDir + "/fimo.txt";
		FileReader r = new FileReader(fimoTxt);
		BufferedReader b = new BufferedReader(r);
		int numLines = 0;
		while(b.ready()) {
			@SuppressWarnings("unused")
			String line = b.readLine();
			numLines++;
		}
		r.close();
		b.close();
		return numLines > 1;
	}
	
	/**
	 * Run DREME on each gene individually
	 * @param dremeExecutable DREME executable
	 * @param batchedMotifPipelineJar BatchedMotifPipeline jar file
	 * @return List of LSF job IDs
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws DrmaaException 
	 */
	private ArrayList<Job> submitDremeOnEachGene(String dremeExecutable, String batchedMotifPipelineJar, Scheduler scheduler) throws IOException, InterruptedException, DrmaaException {
		logger.info("");
		logger.info("Running DREME on each gene individually...");
		ArrayList<Job> jobs = new ArrayList<Job>();
		for(String chr : genes.keySet()) {
			for(Gene gene : genes.get(chr)) {
				makeIndividualGeneBedFiles(gene);
				if(!fileHasData(getIndividualGeneOverlappingFeaturesBedFileName(gene))) {
					logger.warn("Gene " + gene.getName() + " has no overlapping features. Skipping.");
					continue;
				}
				Job job = submitDremeOnGene(gene, dremeExecutable, batchedMotifPipelineJar, scheduler, session);
				jobs.add(job);
			}
		}
		logger.info("All jobs submitted.");
		return jobs;
	}
	
	/**
	 * Run Fimo on each gene individually
	 * @param fimoExecutable Fimo executable
	 * @param batchedMotifPipelineJar BatchedMotifPipeline jar file
	 * @return List of LSF job IDs
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws DrmaaException 
	 */
	private ArrayList<Job> submitFimoOnEachGene(String fimoExecutable, String batchedMotifPipelineJar, Scheduler scheduler) throws IOException, InterruptedException, DrmaaException {
		logger.info("");
		logger.info("Running FIMO on each gene individually...");
		ArrayList<Job> jobs = new ArrayList<Job>();
		for(String chr : genes.keySet()) {
			for(Gene gene : genes.get(chr)) {
				if(!dremeFoundMotif(getJobDescription(gene))) {
					logger.warn("Gene " + gene.getName() + " has no motifs. Skipping.");
					continue;
				}
				Job job = submitFimoOnGene(gene, fimoExecutable, ESTIMATE_FIMO_ALPHA_INDIVIDUAL_GENE, DEFAULT_FIMO_QVAL_THRESHOLD, batchedMotifPipelineJar, scheduler, session);
				jobs.add(job);
			}
		}
		logger.info("All jobs submitted.");
		return jobs;
	}
	
	
	
	/**
	 * Submit job to convert fimo.txt file to bed file for each gene
	 * @param fimo2BedJar Fimo2Bed jar file
	 * @return List of LSF job IDs
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws DrmaaException 
	 */
	private ArrayList<Job> submitFimo2BedOnEachGene(String fimo2BedJar, Scheduler scheduler) throws IOException, InterruptedException, DrmaaException {
		logger.info("");
		logger.info("Creating bed files of FIMO results for each gene individually...");
		ArrayList<Job> jobs = new ArrayList<Job>();
		for(String chr : genes.keySet()) {
			for(Gene gene : genes.get(chr)) {
				String description = getJobDescription(gene);
				if(!dremeFoundMotif(description)) {
					continue;
				}
				if(!fimoFoundMotif(description)) {
					logger.warn("Gene " + gene.getName() + " has no motif occurrences from FIMO. Skipping.");
					continue;
				}
				String bed = getIndividualGeneOverlappingCollapsedFeaturesBedFileName(gene);
				String fimoDir = getFimoDirectory(description);
				String fimoFile = fimoDir + "/fimo.txt";
				String outDir = FIMO_DIR;
				String id = "fimo_" + description;
				Fimo2BedJob f = new Fimo2BedJob(fimoFile, bed, id, outDir, fimo2BedJar, scheduler);
				Job job = f.submitJob(session);
				jobs.add(job);
			}
		}
		return jobs;
	}
	
	private ArrayList<Job> submitFimo2Bed(String fimo2BedJar, Scheduler scheduler) throws IOException, InterruptedException, DrmaaException {
		logger.info("");
		logger.info("Creating bed files of FIMO results for all genes together...");
		ArrayList<Job> jobs = new ArrayList<Job>();
		for(String description : fimoJobDescriptions) {
			if(!fimoFoundMotif(description)) {
				logger.warn("There were no motif occurrences from FIMO. Skipping.");
				return jobs;
			}
			String fimoDir = getFimoDirectory(description);
			String fimoFile = fimoDir + "/fimo.txt";
			String outDir = FIMO_DIR;
			String id = description;
			Fimo2BedJob f = new Fimo2BedJob(fimoFile, geneBed, id, outDir, fimo2BedJar, scheduler);
			Job job = f.submitJob(session);
			jobs.add(job);
		}
		return jobs;
	}
	
	
	private static boolean fileHasData(String fileName) throws IOException {
		FileReader r = new FileReader(fileName);
		BufferedReader b = new BufferedReader(r);
		boolean hasData = false;
		while(b.ready()) {
			String line = b.readLine();
			if(line.length() > 0) {
				hasData = true;
			}
		}
		r.close();
		b.close();
		return hasData;
	}
	
	private static String getDremeDirectory(String jobDescription) {
		return DREME_DIR + "/" + jobDescription;
	}
	
	private static String getFimoDirectory(String jobDescription) {
		return FIMO_DIR + "/" + jobDescription;
	}
	
	/**
	 * Run DREME on features
	 * @param dremeExecutable DREME executable
	 * @param additionalDremeOptions Additional DREME options
	 * @param jobDescription Job description
	 * @param lsfQueue LSF queue
	 * @param memoryRequest Memory request in Gb
	 * @return LSF job ID
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws DrmaaException 
	 */
	protected Job submitDreme(String dremeExecutable, String additionalDremeOptions, String jobDescription, String lsfQueue, int memoryRequest, Scheduler scheduler) throws IOException, InterruptedException, DrmaaException {
		logger.info("");
		logger.info("Running DREME...");
		DremeJob d = new DremeJob(dremeExecutable, collapsedFeatureFasta, genesMinusCollapsedFeaturesWindowFasta, jobDescription, getDremeDirectory(jobDescription), additionalDremeOptions, scheduler);
		Job job = d.submitJob(lsfQueue, memoryRequest, session);
		dremeJobDescriptions.add(jobDescription);
		return job;
	}
	
	/**
	 * Get motif xml file generated by Dreme
	 * @param jobDescription Job description
	 * @return File path
	 */
	protected static String getMotifFile(String jobDescription) {
		return getDremeDirectory(jobDescription) + "/dreme.xml";
	}
	
	/**
	 * Run FIMO
	 * @param fimoExecutable Fimo executable
	 * @param fimoOptionAlpha Expected proportion of sequences that harbor the motif
	 * @param fimoOptionQvalThresh Q value threshold
	 * @param additionalFimoOptions Additional Fimo options
	 * @param jobDescription Job description
	 * @param lsfQueue LSF queue
	 * @param memoryRequest Memory request in Gb
	 * @return LSF job ID
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws DrmaaException 
	 */
	protected Job submitFimo(String fimoExecutable, double fimoOptionAlpha, double fimoOptionQvalThresh, String additionalFimoOptions, String jobDescription, String lsfQueue, int memoryRequest, Scheduler scheduler) throws IOException, InterruptedException, DrmaaException {
		logger.info("");
		logger.info("Running FIMO...");
		FimoJob f = new FimoJob(fimoExecutable, collapsedFeatureFasta, getMotifFile(jobDescription), getFimoDirectory(jobDescription), fimoOptionAlpha, fimoOptionQvalThresh, additionalFimoOptions, jobDescription);
		Job job = f.submitJob(lsfQueue, memoryRequest, scheduler, session);
		fimoJobDescriptions.add(jobDescription);
		return job;
	}
	
	/**
	 * Run FIMO on all features with motifs already found in individual genes
	 * @param fimoExecutable Fimo executable
	 * @param fimoOptionAlpha Expected proportion of sequences that harbor the motif
	 * @param fimoOptionQvalThresh Q value threshold
	 * @param additionalFimoOptions Additional Fimo options
	 * @return List of LSF job IDs
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws DrmaaException 
	 */
	private ArrayList<Job> submitFimoOnAllGenesWithIndividualGeneMotifs(String fimoExecutable, double fimoOptionAlpha, double fimoOptionQvalThresh, String additionalFimoOptions, Scheduler scheduler) throws IOException, InterruptedException, DrmaaException {
		logger.info("");
		logger.info("Using FIMO to search for individual gene motifs among all genes...");
		ArrayList<Job> jobs = new ArrayList<Job>();
		for(String chr : genes.keySet()) {
			for(Gene gene : genes.get(chr)) {
				String geneJobDescription = getJobDescription(gene);
				if(!dremeFoundMotif(geneJobDescription)) {
					continue;
				}
				logger.info("Searching for motifs from gene " + gene.getName());
				String jobDescription = ALL_FEATURES_JOB_DESCRIPTION + "_" + gene.getName() + "_motif";
				FimoJob f = new FimoJob(fimoExecutable, collapsedFeatureFasta, getMotifFile(geneJobDescription), getFimoDirectory(jobDescription), fimoOptionAlpha, fimoOptionQvalThresh, additionalFimoOptions, jobDescription);
				Job job = f.submitJob("week", 4, scheduler, session);
				fimoJobDescriptions.add(jobDescription);
				jobs.add(job);
			}
		}
		logger.info("All jobs submitted.");
		return jobs;
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 * @throws InterruptedException 
	 * @throws DrmaaException 
	 */
	public static void main(String[] args) throws IOException, InterruptedException, DrmaaException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-gb", "Gene bed file", true);
		p.addStringArg("-fb", "Feature bed file", true);
		p.addStringArg("-g", "Genome fasta file", true);
		p.addStringArg("-d", "Dreme executable", true);
		p.addStringArg("-bj", "BatchedMotifPipeline jar", true);
		p.addStringArg("-f", "Fimo executable", true);
		p.addStringArg("-f2b", "Fimo2Bed jar", true);
		p.addStringArg("-s", "Name of scheduler e.g. LSF or OGS", true);
		p.parse(args);
		String geneBed = p.getStringArg("-gb");
		String featureBed = p.getStringArg("-fb");
		String genomeFasta = p.getStringArg("-g");
		String dremeExecutable = p.getStringArg("-d");
		String batchedJar = p.getStringArg("-bj");
		String fimoExecutable = p.getStringArg("-f");
		String fimo2bedJar = p.getStringArg("-f2b");
		Scheduler scheduler = Scheduler.fromString(p.getStringArg("-s"));
		Session drmaaSession = scheduler.equals(Scheduler.OGS) ? OGSUtils.getDrmaaSession() : null;
		
		MotifPipeline m = new MotifPipeline(geneBed, featureBed, genomeFasta, drmaaSession);

		// Run DREME on all genes together
		ArrayList<Job> allGenesDremeJob = new ArrayList<Job>();
		Job job = m.submitDreme(dremeExecutable, DremeJob.ADDITIONAL_OPTIONS, ALL_FEATURES_JOB_DESCRIPTION, "week", 8, scheduler);
		allGenesDremeJob.add(job);

		// Run DREME on each gene individually
		ArrayList<Job> indGeneDremeJobIDs = m.submitDremeOnEachGene(dremeExecutable, batchedJar, scheduler);
		// Wait for individual gene dreme jobs to finish
		logger.info("");
		logger.info("Waiting for individual gene DREME jobs to finish...");
		JobUtils.waitForAll(indGeneDremeJobIDs);
		
		// Run FIMO on all genes with motifs found in individual genes
		ArrayList<Job> allGenesFimoJobs = new ArrayList<Job>();
		ArrayList<Job> fJobs = m.submitFimoOnAllGenesWithIndividualGeneMotifs(fimoExecutable, 0.001, DEFAULT_FIMO_QVAL_THRESHOLD, FimoJob.ADDITIONAL_OPTIONS, scheduler);
		allGenesFimoJobs.addAll(fJobs);
		
		// Run FIMO on each gene individually
		ArrayList<Job> indGeneFimoJobs = m.submitFimoOnEachGene(fimoExecutable, batchedJar, scheduler);
		// Wait for individual gene fimo jobs to finish
		logger.info("");
		logger.info("Waiting for individual gene FIMO jobs to finish...");
		JobUtils.waitForAll(indGeneFimoJobs);

		// Run Fimo2Bed on individual genes
		logger.info("");
		logger.info("Converting FIMO output to bed tracks for individual genes...");
		ArrayList<Job> indGeneFimo2BedJobs = m.submitFimo2BedOnEachGene(fimo2bedJar, scheduler);
		
		// Wait for individual gene Fimo2Bed jobs to finish
		logger.info("");
		logger.info("Waiting for individual gene Fimo2Bed jobs to finish...");
		JobUtils.waitForAll(indGeneFimo2BedJobs);
		
		// Wait for all genes DREME job to finish
		logger.info("");
		logger.info("Waiting for DREME to finish on all genes...");
		JobUtils.waitForAll(allGenesDremeJob);
		
		// Run FIMO on all genes
		Job fJob = m.submitFimo(fimoExecutable, 0.001, DEFAULT_FIMO_QVAL_THRESHOLD, FimoJob.ADDITIONAL_OPTIONS, ALL_FEATURES_JOB_DESCRIPTION, "week", 8, scheduler);
		allGenesFimoJobs.add(fJob);		
		
		// Waiting for all genes FIMO jobs to finish
		logger.info("");
		logger.info("Waiting for FIMO jobs to finish on all genes...");
		JobUtils.waitForAll(allGenesFimoJobs);
		
		// Run Fimo2Bed on all genes
		logger.info("");
		logger.info("Converting FIMO output to bed tracks for all genes...");
		ArrayList<Job> allGeneFimo2BedJobs = m.submitFimo2Bed(fimo2bedJar, scheduler);
		
		// Waiting for all genes Fimo2Bed job to finish
		logger.info("");
		logger.info("Waiting for Fimo2Bed to finish on all genes...");
		JobUtils.waitForAll(allGeneFimo2BedJobs);
		
		logger.info("");
		logger.info("All done.");

	}

}
