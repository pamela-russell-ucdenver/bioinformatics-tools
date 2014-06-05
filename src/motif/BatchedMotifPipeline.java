/**
 * 
 */
package motif;

import java.io.IOException;
import java.util.ArrayList;

import org.ggf.drmaa.DrmaaException;

import nextgen.core.job.Job;
import nextgen.core.job.JobUtils;
import nextgen.core.pipeline.Scheduler;


import broad.core.parser.CommandLineParser;

/**
 * @author prussell
 *
 */
public class BatchedMotifPipeline extends MotifPipeline {

	// test comment

	/**
	 * @param geneBedFile
	 * @param featureBedFile
	 * @param genomeFastaFile
	 * @throws IOException
	 */
	public BatchedMotifPipeline(String geneBedFile, String featureBedFile, String genomeFastaFile) throws IOException {
		super(geneBedFile, featureBedFile, genomeFastaFile);
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
		p.addStringArg("-d", "Dreme executable", false, null);
		p.addStringArg("-f", "Fimo executable", false, null);
		p.addStringArg("-de", "Description", true);
		p.addDoubleArg("-a", "Fimo alpha option: expected proportion of sequences with the motif", true);
		p.addDoubleArg("-q", "Q value threshold", true);
		p.addStringArg("-s", "Scheduler e.g. LSF or OGS", true);
		p.parse(args);
		String geneBed = p.getStringArg("-gb");
		String featureBed = p.getStringArg("-fb");
		String genomeFasta = p.getStringArg("-g");
		String dremeExecutable = p.getStringArg("-d");
		String description = p.getStringArg("-de");
		String fimoExecutable = p.getStringArg("-f");
		double alpha = p.getDoubleArg("-a");
		double qvalThresh = p.getDoubleArg("-q");
		Scheduler scheduler = Scheduler.fromString(p.getStringArg("-s"));

		BatchedMotifPipeline m = new BatchedMotifPipeline(geneBed, featureBed, genomeFasta);

		if(dremeExecutable != null) {
			ArrayList<Job> jobs = new ArrayList<Job>();
			// Run DREME
			Job job = m.submitDreme(dremeExecutable, DremeJob.ADDITIONAL_OPTIONS, description, "hour", 8, scheduler);
			jobs.add(job);
			// Wait for DREME job to finish
			JobUtils.waitForAll(jobs);
		}

		if(fimoExecutable != null) {
			ArrayList<Job> jobs = new ArrayList<Job>();
			// Run Fimo
			Job job = m.submitFimo(fimoExecutable, alpha, qvalThresh, FimoJob.ADDITIONAL_OPTIONS, description, "hour", 4, scheduler);
			jobs.add(job);
			// Wait for Fimo job to finish
			JobUtils.waitForAll(jobs);
		}
		
	}

}
