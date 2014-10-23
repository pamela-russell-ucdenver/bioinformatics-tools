/**
 * 
 */
package motif;

import java.io.IOException;

import org.ggf.drmaa.DrmaaException;
import org.ggf.drmaa.Session;

import guttmanlab.core.pipeline.Job;
import guttmanlab.core.pipeline.LSFJob;
import guttmanlab.core.pipeline.OGSJob;
import guttmanlab.core.pipeline.Scheduler;



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
	 * @param drmaaSession Active DRMAA session or null if not using OGS
	 * @return Job
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws DrmaaException 
	 */
	public Job submitJob(Session drmaaSession) throws IOException, InterruptedException, DrmaaException {
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
		case OGS:
			OGSJob ogsjob = new OGSJob(drmaaSession, cmmd);
			ogsjob.submit();
			return ogsjob;
		default:
			throw new IllegalArgumentException("Scheduler " + scheduler.toString() + " not supported.");
		}
	}
	
}
