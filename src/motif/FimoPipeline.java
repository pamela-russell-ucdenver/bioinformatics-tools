/**
 * 
 */
package motif;

import general.CommandLineParser;

import java.io.IOException;



import nextgen.core.pipeline.util.OGSUtils;

import org.apache.log4j.Logger;
import org.ggf.drmaa.DrmaaException;
import org.ggf.drmaa.Session;

import pipeline.Job;
import pipeline.Scheduler;


/**
 * @author prussell
 *
 */
public class FimoPipeline {

	private String jobDescription;
	private String outDir;
	private String fimoTxtFile;
	private static Logger logger = Logger.getLogger(FimoPipeline.class.getName());
	private Session session;
	
	
	private FimoPipeline(String description, String outDirectory, Session drmaaSession) {
		jobDescription = description;
		outDir = outDirectory;
		fimoTxtFile = outDir + "/fimo.txt";
		session = drmaaSession;
	}
	
	private void runFimo(String fimoExecutable, String sequenceFasta, String motifXml, double fimoOptionAlpha, double fimoOptionQvalueThreshold, String addlOptions, String queue, int memoryRequest, Scheduler scheduler) throws IOException, InterruptedException, DrmaaException {
		logger.info("");
		logger.info("Running FIMO...");
		FimoJob fj = new FimoJob(fimoExecutable, sequenceFasta, motifXml, outDir, fimoOptionAlpha, fimoOptionQvalueThreshold, addlOptions, jobDescription);
		Job job = fj.submitJob(queue, memoryRequest, scheduler, session);
		logger.info("Waiting for FIMO job to finish...");
		job.waitFor();
		logger.info("Done running FIMO.");
	}
	
	private void runFimo2Bed(String bedAnnotation, String fimo2BedJar, Scheduler scheduler) throws IOException, InterruptedException, DrmaaException {
		logger.info("");
		logger.info("Running Fimo2Bed...");
		Fimo2BedJob f2b = new Fimo2BedJob(fimoTxtFile, bedAnnotation, jobDescription, outDir, fimo2BedJar, scheduler);
		Job job = f2b.submitJob(session);
		logger.info("Waiting for Fimo2Bed job to finish...");
		job.waitFor();
		logger.info("Done running Fimo2Bed.");
	}
	
	/**
	 * @param args
	 * @throws InterruptedException 
	 * @throws IOException 
	 * @throws DrmaaException 
	 */
	public static void main(String[] args) throws IOException, InterruptedException, DrmaaException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-fe", "Fimo executable file", true);
		p.addStringArg("-sf", "Sequence fasta file", true);
		p.addStringArg("-mx", "Motif xml file", true);
		p.addDoubleArg("-fa", "Fimo alpha option", true);
		p.addDoubleArg("-fq", "Fimo Q value threshold", true);
		p.addStringArg("-q", "LSF queue", false, "hour");
		p.addIntArg("-mr", "Memory request", false, 4);
		p.addStringArg("-ba", "Bed annotation", true);
		p.addStringArg("-f2b", "Fimo2Bed jar file", true);
		p.addStringArg("-d", "Description", true);
		p.addStringArg("-o", "Output directory", true);
		p.addStringArg("-s", "Scheduler e.g. LSF or OGS", true);
		p.parse(args);
		String fimoExecutable = p.getStringArg("-fe");
		String sequenceFasta = p.getStringArg("-sf");
		String motifXml = p.getStringArg("-mx");
		double fimoOptionAlpha = p.getDoubleArg("-fa");
		double fimoOptionQvalueThreshold = p.getDoubleArg("-fq");
		String queue = p.getStringArg("-q");
		int memoryRequest = p.getIntArg("-mr");
		String bedAnnotation = p.getStringArg("-ba");
		String fimo2BedJar = p.getStringArg("-f2b");
		String description = p.getStringArg("-d");
		String outDirectory = p.getStringArg("-o");
		Scheduler scheduler = Scheduler.fromString(p.getStringArg("-s"));
		
		Session drmaaSession = scheduler.equals(Scheduler.OGS) ? OGSUtils.getDrmaaSession() : null;
		
		FimoPipeline fp = new FimoPipeline(description, outDirectory, drmaaSession);
		fp.runFimo(fimoExecutable, sequenceFasta, motifXml, fimoOptionAlpha, fimoOptionQvalueThreshold, FimoJob.ADDITIONAL_OPTIONS, queue, memoryRequest, scheduler);
		fp.runFimo2Bed(bedAnnotation, fimo2BedJar, scheduler);
		
		logger.info("");
		logger.info("All done.");
		
	}

}
