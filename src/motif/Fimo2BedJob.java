/**
 * 
 */
package motif;

import java.io.IOException;

import nextgen.core.job.Job;
import nextgen.core.job.LSFJob;
import nextgen.core.pipeline.Scheduler;


/**
 * @author prussell
 *
 */
public class Fimo2BedJob {

	private String bedAnnotation;
	private String fimoFile;
	private String outDir;
	private String id;
	private String fimo2bedJar;
	private Scheduler scheduler;
	
	/**
	 * @param fimoTxtFile
	 * @param bedFile
	 * @param description
	 * @param outputDir
	 * @param fimo2BedJar
	 * @param sched
	 */
	public Fimo2BedJob(String fimoTxtFile, String bedFile, String description, String outputDir, String fimo2BedJar, Scheduler sched) {
		bedAnnotation = bedFile;
		fimoFile = fimoTxtFile;
		outDir = outputDir;
		id = description;
		fimo2bedJar = fimo2BedJar;
		scheduler = sched;
	}
	
	/**
	 * @return Job
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public Job submitJob() throws IOException, InterruptedException {
		String cmmd = "java -jar -Xmx3g -Xms2g -Xmn1g " + fimo2bedJar + " ";
		cmmd += "-b " + bedAnnotation + " ";
		cmmd += "-f " + fimoFile + " ";
		cmmd += "-id " + id + " ";
		cmmd += "-o " + outDir + " ";
		switch(scheduler) {
		case LSF:
			String lsfJobID = id + "_" + System.currentTimeMillis();
			String bsubOut = outDir + "/fimo_to_bed_" + lsfJobID + ".bsub";
			LSFJob job = new LSFJob(Runtime.getRuntime(), lsfJobID, cmmd, bsubOut, "hour", 4);
			job.submit();
			return job;
		default:
			throw new IllegalArgumentException("Scheduler " + scheduler.toString() + " not supported.");
		}
	}
	
}
