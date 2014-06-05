/**
 * 
 */
package tests;

import java.io.IOException;
import java.util.ArrayList;

import nextgen.core.job.Job;
import nextgen.core.job.JobUtils;
import nextgen.core.job.OGSJob;
import nextgen.core.pipeline.util.OGSUtils;

import org.ggf.drmaa.DrmaaException;
import org.ggf.drmaa.Session;

/**
 * @author prussell
 *
 */
public class TestOGSJob {
	
	private Session drmaaSession;
	
	private TestOGSJob(Session session) {
		drmaaSession = session;
	}
	
	private OGSJob runCommand(String cmmd, String jobName, boolean deleteScriptFile, boolean waitFor) throws DrmaaException, IOException, InterruptedException {
		OGSJob job = new OGSJob(drmaaSession, cmmd, deleteScriptFile, jobName, null);
		job.submit();
		if(waitFor) job.waitFor();
		return job;
	}
	
	
	/**
	 * @param args
	 * @throws DrmaaException 
	 * @throws InterruptedException 
	 * @throws IOException 
	 */
	public static void main(String[] args) throws DrmaaException, IOException, InterruptedException {
		
		TestOGSJob t = new TestOGSJob(OGSUtils.getDrmaaSession());
		
		OGSJob pwdJob = t.runCommand("pwd", "pwd", true, false);
		OGSJob lsJob = t.runCommand("ls", "ls", true, false);
		String cmmd = "/storage/Software/Fastx/fastx_0.0.13/fastx_clipper -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -Q 33 -n -i 4K_FA_2.fq -o 4K_FA_2.fq.clipped_tmp";
		//OGSJob fastxJob = t.runCommand(cmmd, "fastx_from_string", false, false);
		
		OGSJob ogsJob = new OGSJob(t.drmaaSession, cmmd, true, "fastx_clip_adapters", null);
		ogsJob.submit();
		ArrayList<Job> jobs = new ArrayList<Job>();

		jobs.add(pwdJob);
		jobs.add(lsJob);
		jobs.add(ogsJob);
		//jobs.add(fastxJob);
		JobUtils.waitForAll(jobs);
		
	}

}
