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
public class DremeJob {

	/*

    /seq/lincRNA/Pam/Software/meme_4.9.0/scripts/dreme [options]

    -o  <directory>         create the specified output directory 
                            and write all output to files in that directory
    -oc <directory>         create the specified output directory 
                            overwritting it if it already exists;
                            default: create dreme_out in the currrent
                            working directory
    -p <filename>           positive sequence file name (required)
    -n <filename>           negative sequence file name (optional);
                            default: the positive sequences are shuffled
                            to create the negative set if -n is not used
    -norc                   search given strand only for motifs (not reverse complement)
    -e <ethresh>            stop if motif E-value > <ethresh>;
                            default: 0.05
    -m <m>                  stop if <m> motifs have been output;
                            default: only stop at E-value threshold
    -t <seconds>            stop if the specified time has elapsed;
                            default: only stop at E-value threshold
    -g <ngen>               number of REs to generalize; default: 100
                            Hint: Increasing <ngen> will make the motif
                            search more thoroughly at some cost in speed.
    -s <seed>               seed for shuffling sequences; ignored
                            if -n <filename> given; default: 1
    -v <verbosity>          1..5 for varying degrees of extra output
                            default: 2
    -png                    create PNG logos
    -eps                    create EPS logos
    -desc <description>     store the description in the output;
                            default: no description
    -dfile <filename>       acts like -desc but reads the description from
                            the specified file; allows characters that would 
                            otherwise have to be escaped; 
                            default: no description
    -h                      print this usage message

-----------------------Setting Core Motif Width---------------------------------
                   Hint: The defaults are pretty good; making k larger
                         than 8 slows DREME down with little other effect.
                         Use these if you just want motifs shorter than 8.
--------------------------------------------------------------------------------
    -mink <mink>            minimum width of core motif; default 3
    -maxk <maxk>            maximum width of core motif; default 8
    -k <k>                  sets mink=maxk=<k>

	 */

	private String executable;
	private String positiveSeqFasta;
	private String negativeSeqFasta;
	private String outDir;
	private String description;
	private String addlOptions;
	private String lsfJobID;
	private Scheduler scheduler;

	/**
	 * Additional DREME options
	 */
	public static String ADDITIONAL_OPTIONS = "-norc -e 0.05 -eps -png -mink 3 -maxk 20 -v 5 -g 1000";
	
	/**
	 * @param executableFileName Dreme executable python script
	 * @param positiveSequenceFasta Fasta file of positive sequences
	 * @param negativeSequenceFasta Fasta file of negative sequences
	 * @param outDirectory Output directory
	 * @param jobDescription Description
	 * @param sched Scheduler
	 */
	public DremeJob(String executableFileName, String positiveSequenceFasta, String negativeSequenceFasta, String jobDescription, String outDirectory, Scheduler sched) {
		this(executableFileName, positiveSequenceFasta, negativeSequenceFasta, jobDescription, outDirectory, "", sched);
	}
	
	/**
	 * @param executableFileName Dreme executable python script
	 * @param positiveSequenceFasta Fasta file of positive sequences
	 * @param negativeSequenceFasta Fasta file of negative sequences
	 * @param jobDescription Description
	 * @param outDirectory Output directory
	 * @param additionalDremeOptions Additional dreme flags and values
	 * @param sched Scheduler
	 */
	public DremeJob(String executableFileName, String positiveSequenceFasta, String negativeSequenceFasta, String jobDescription, String outDirectory, String additionalDremeOptions, Scheduler sched) {
		executable = executableFileName;
		positiveSeqFasta = positiveSequenceFasta;
		negativeSeqFasta = negativeSequenceFasta;
		description = jobDescription;
		outDir = outDirectory;
		addlOptions = additionalDremeOptions;
		scheduler = sched;
	}
		
	/**
	 * Submit the job to LSF
	 * @param queue Queue name
	 * @param memoryRequest Memory request in gigabytes
	 * @param drmaaSession Active DRMAA session or null if not using OGS
	 * @return Job ID
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws DrmaaException 
	 */
	public Job submitJob(String queue, int memoryRequest, Session drmaaSession) throws IOException, InterruptedException, DrmaaException {
		lsfJobID = description + "_" + System.currentTimeMillis();
		String bsubOut = outDir + "/" + lsfJobID + ".bsub";
		String command = "python " + executable + " ";
		command += "-oc " + outDir + " ";
		command += "-p " + positiveSeqFasta + " ";
		command += "-n " + negativeSeqFasta + " ";
		command += "-desc " + description + " ";
		command += addlOptions;
		
		switch(scheduler) {
		case LSF:
			LSFJob job = new LSFJob(Runtime.getRuntime(), lsfJobID, command, bsubOut, queue, memoryRequest);
			job.submit();
			return job;
		case OGS:
			OGSJob ogsjob = new OGSJob(drmaaSession, command);
			ogsjob.submit();
			return ogsjob;
		default:
			throw new IllegalArgumentException("Scheduler " + scheduler + " not supported.");
		}
	}
	
	
}
