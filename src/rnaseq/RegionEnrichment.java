/**
 * 
 */
package rnaseq;

import guttmanlab.core.util.CommandLineParser;

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

import guttmanlab.core.pipeline.Job;
import guttmanlab.core.pipeline.JobUtils;
import guttmanlab.core.pipeline.LSFJob;
import guttmanlab.core.pipeline.Scheduler;

import broad.pda.annotation.BEDFileParser;

import nextgen.core.annotation.Gene;
import nextgen.core.coordinatesystem.TranscriptomeSpace;
import nextgen.core.model.AlignmentModel;
import nextgen.core.normalize.CrossSampleTranscriptAverageNormalization;
import nextgen.core.normalize.TranscriptAverageNormalization;

/**
 * @author prussell
 *
 */
public class RegionEnrichment {
	
	TranscriptAverageNormalization enrichments;
	CrossSampleTranscriptAverageNormalization comparativeEnrichments;
	private Map<String, Collection<Gene>> genes;
	private static Logger logger = Logger.getLogger(RegionEnrichment.class.getName());
	private Collection<String> chrs;
	private Map<String, Collection<Gene>> features;
	
	private RegionEnrichment(String bamFileDataToNormalize, String bamFileDataToNormalizeTo, String geneBedFile) throws IOException {
		this(bamFileDataToNormalize, bamFileDataToNormalizeTo, geneBedFile, null);
	}
	
	private RegionEnrichment(String bamFileDataToNormalize, String bamFileDataToNormalizeTo, String geneBedFile, String chr) throws IOException {
		this(bamFileDataToNormalize, bamFileDataToNormalizeTo, geneBedFile, null, chr);
	}
	
	private RegionEnrichment(String bamFileDataToNormalize, String bamFileDataToNormalizeTo, String geneBedFile, String regionBedFile, String chr) throws IOException {
		logger.info("");
		logger.info("First normalizing counts within transcript, then normalizing " + bamFileDataToNormalize + " to " + bamFileDataToNormalizeTo + ".");
		genes = BEDFileParser.loadDataByChr(new File(geneBedFile));
		if(regionBedFile != null) features = BEDFileParser.loadDataByChr(new File(regionBedFile));
		chrs = new ArrayList<String>();
		if(chr != null) {
			chrs.add(chr);
		} else {
			chrs.addAll(genes.keySet());
		}
		logger.info("");
		logger.info("Using chromosomes:");
		for(String c : chrs) {
			logger.info(c);
		}
		TranscriptomeSpace t = new TranscriptomeSpace(genes);
		AlignmentModel dataToNormalize = new AlignmentModel(bamFileDataToNormalize, t);
		AlignmentModel dataToNormalizeTo = new AlignmentModel(bamFileDataToNormalizeTo, t);
		enrichments = regionBedFile == null ? new TranscriptAverageNormalization(dataToNormalize) : new TranscriptAverageNormalization(dataToNormalize, regionBedFile, geneBedFile);
		TranscriptAverageNormalization normalizedDataToNormalizeTo = regionBedFile == null ? new TranscriptAverageNormalization(dataToNormalizeTo) : new TranscriptAverageNormalization(dataToNormalizeTo, regionBedFile, geneBedFile);
		comparativeEnrichments = new CrossSampleTranscriptAverageNormalization(enrichments, normalizedDataToNormalizeTo);
	}
	
	
	
	private static void batchRunAllChromosomes(CommandLineParser p) throws IOException, InterruptedException, DrmaaException {
		logger.info("");
		logger.info("Running each chromosome on LSF...");
		String jar = p.getStringArg("-j");
		if(jar == null) {
			throw new IllegalArgumentException("Must specify jar file.");
		}
		String bam1 = p.getStringArg("-b1");
		String bam2 = p.getStringArg("-b2");
		String bed = p.getStringArg("-g");
		String wigPrefix = p.getStringArg("-w");
		String wigToBigWig = p.getStringArg("-u");
		String sizeFile = p.getStringArg("-s");
		String featureBed = p.getStringArg("-f");
		String outFeatureBedPrefix = p.getStringArg("-fb");
		Scheduler scheduler = Scheduler.fromString(p.getStringArg("-sc"));
		RegionEnrichment c = new RegionEnrichment(bam1, bam2, bed, featureBed, null);
		Collection<String> chrList = c.chrs;
		Collection<String> argsToIgnore = new ArrayList<String>();
		argsToIgnore.add("-w");
		argsToIgnore.add("-c");
		argsToIgnore.add("-bc");
		argsToIgnore.add("-j");
		argsToIgnore.add("-s");
		argsToIgnore.add("-u");
		argsToIgnore.add("-fb");
		String argString = p.getArgString(argsToIgnore);
		ArrayList<Job> jobs = new ArrayList<Job>();
		ArrayList<String> chrWigs = new ArrayList<String>();
		ArrayList<String> chrFeatureBeds = new ArrayList<String>();
		for(String chr : chrList) {
			String cmmd =  "java -jar -Xmx30g -Xms25g -Xmn20g " + jar + " " + argString + " -c " + chr + " ";
			if(wigPrefix != null) {
				String outWig = wigPrefix + "_" + chr;
				cmmd += "-w " + outWig + " ";
				chrWigs.add(outWig + ".wig");
			}
			if(outFeatureBedPrefix != null) {
				if(featureBed == null) {
					throw new IllegalArgumentException("Must provide bed file of features.");
				}
				String outBed = outFeatureBedPrefix + "_" + chr + ".bed";
				cmmd += " -fb " + outBed + " ";
				chrFeatureBeds.add(outBed);
			}
			switch(scheduler) {
			case LSF:
				String jobID = Long.valueOf(System.currentTimeMillis()).toString();
				String bsubOut = chr + "_" + jobID + ".bsub";
				logger.info("Running command " + cmmd);
				LSFJob job = new LSFJob(Runtime.getRuntime(), jobID, cmmd, bsubOut, "week", 32);
				job.submit();
				jobs.add(job);
				break;
			default:
				throw new IllegalArgumentException("Scheduler " + scheduler.toString() + " not supported.");
			}
		}
		logger.info("");
		logger.info("Waiting for jobs to finish...");
		JobUtils.waitForAll(jobs);
		logger.info("");
		logger.info("All jobs done.");
		if(outFeatureBedPrefix != null) {
			String outFullBed = outFeatureBedPrefix + ".bed";
			logger.info("");
			logger.info("Writing combined bed file to " + outFullBed + "...");
			FileWriter w = new FileWriter(outFullBed);
			for(String chrBed : chrFeatureBeds) {
				logger.info(chrBed);
				FileReader r = new FileReader(chrBed);
				BufferedReader b = new BufferedReader(r);
				while(b.ready()) {
					w.write(b.readLine() + "\n");
				}
				r.close();
				b.close();
			}
			w.close();
			logger.info("Done writing combined bed file.");			
		}
		if(wigPrefix != null) {
			String combinedWig = wigPrefix + ".wig";
			logger.info("");
			logger.info("Writing combined wig file to " + combinedWig + "...");
			FileWriter w = new FileWriter(combinedWig);
			for(String chrWig : chrWigs) {
				logger.info(chrWig);
				FileReader r = new FileReader(chrWig);
				BufferedReader b = new BufferedReader(r);
				while(b.ready()) {
					w.write(b.readLine() + "\n");
				}
				r.close();
				b.close();
			}
			w.close();
			logger.info("Done writing combined wig file.");
			if(wigToBigWig != null && sizeFile != null) {
				NormalizationUtils.makeBigWig(wigPrefix, wigToBigWig, sizeFile);
			}
		}

	}
	
	/**
	 * @param args
	 * @throws IOException 
	 * @throws InterruptedException 
	 * @throws DrmaaException 
	 */
	public static void main(String[] args) throws IOException, InterruptedException, DrmaaException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-b1", "Numerator bam file", true);
		p.addStringArg("-b2", "Denominator bam file", true);
		p.addStringArg("-g", "Gene bed file", true);
		p.addStringArg("-w", "Output wig file prefix", false, null);
		p.addStringArg("-c", "Single chromosome to write", false, null);
		p.addBooleanArg("-bc", "For wig file, run each chromosome separately on LSF. Requires -w.", false, false);
		p.addStringArg("-j", "Jar file for batched run", false, null);
		p.addStringArg("-s", "Chromosome size file to make bigwig", false, null);
		p.addStringArg("-u", "UCSC wigToBigWig executable", false, null);
		p.addStringArg("-f", "Input feature bed file", false, null);
		p.addStringArg("-cfb", "Prefix for output bed file of comparative feature enrichments", false, null);
		p.addStringArg("-gfb", "Prefix for output bed file of feature enrichments over gene", false, null);
		p.addStringArg("-sc", "Scheduler e.g. LSF or OGS", true);
		p.parse(args);
		String bam1 = p.getStringArg("-b1");
		String bam2 = p.getStringArg("-b2");
		String bed = p.getStringArg("-g");
		String wigPrefix = p.getStringArg("-w");
		String chr = p.getStringArg("-c");
		boolean batch = p.getBooleanArg("-bc");
		String chrSizeFile = p.getStringArg("-s");
		String wigToBigWig = p.getStringArg("-u");
		String featureBed = p.getStringArg("-f");
		String outComparisonBed = p.getStringArg("-cfb");
		String outEnrichmentOverGeneBed = p.getStringArg("-gfb");
		
		if(!batch) {
			RegionEnrichment c = new RegionEnrichment(bam1, bam2, bed, featureBed, chr);
			if(outComparisonBed != null) {
				logger.info("");
				logger.info("Writing bed file of comparative feature enrichments...");
				if(featureBed == null) {
					throw new IllegalArgumentException("Must provide bed file of features.");
				}				
				NormalizationUtils.writeFeatureEnrichmentBed(c.features, c.comparativeEnrichments, outComparisonBed);
				if(wigPrefix != null) NormalizationUtils.writePositionLevelDataAllGenes(c.comparativeEnrichments, c.features, wigPrefix, wigToBigWig, chrSizeFile);
			}
			if(outEnrichmentOverGeneBed != null) {
				logger.info("");
				logger.info("Writing bed file of feature enrichments over gene...");
				if(featureBed == null) {
					throw new IllegalArgumentException("Must provide bed file of features.");
				}
				NormalizationUtils.writeFeatureEnrichmentBed(c.features, c.enrichments, outEnrichmentOverGeneBed);
				if(wigPrefix != null) NormalizationUtils.writePositionLevelDataAllGenes(c.enrichments, c.features, wigPrefix, wigToBigWig, chrSizeFile);
			}
		} else {
			batchRunAllChromosomes(p);
		}
		
		logger.info("");
		logger.info("All done.");
		
	}

}
