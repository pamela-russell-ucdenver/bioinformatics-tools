package test;

import guttmanlab.core.util.CommandLineParser;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;

import org.apache.log4j.Logger;
import org.ggf.drmaa.DrmaaException;
import org.ggf.drmaa.Session;

import guttmanlab.core.pipeline.Job;
import guttmanlab.core.pipeline.JobUtils;
import guttmanlab.core.pipeline.OGSJob;

import nextgen.core.pipeline.util.OGSUtils;


public class DrmaaSleepTest {
	
	public static Logger logger = Logger.getLogger(DrmaaSleepTest.class.getName());
	
	/**
	 * @param args
	 * @throws DrmaaException 
	 * @throws IOException 
	 * @throws InterruptedException 
	 */
	public static void main(String[] args) throws DrmaaException, IOException, InterruptedException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-s", "Sleep.jar file", true);
		p.addIntArg("-n", "Number of jobs to submit", true);
		p.addIntArg("-m", "Number of minutes for each job to run", true);
		p.addStringArg("-e", "Email address", true);
		
		p.parse(args);
		String jar = p.getStringArg("-s");
		int numJobs = p.getIntArg("-n");
		int mins = p.getIntArg("-m");
		String email = p.getStringArg("-e");
		
		String cmmd = "java -jar " + jar + " -m " + mins;
		
		Collection<Job> jobs = new ArrayList<Job>();
		Session drmaaSession = OGSUtils.getDrmaaSession();
		
		for(int i = 0; i < numJobs; i++) {
			OGSJob job = new OGSJob(drmaaSession, cmmd, "job" + i, email);
			job.submit();
			jobs.add(job);
		}
		
		JobUtils.waitForAll(jobs);
		
		logger.info("");
		logger.info("All done.");
		
	}

}
