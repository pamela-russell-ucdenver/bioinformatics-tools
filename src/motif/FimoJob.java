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
public class FimoJob {

	/*

The name FIMO stands for "Find Individual Motif Occurences." The program searches a sequence database for occurrences of known motifs, treating each motif independently. The program uses a dynamic programming algorithm to convert log-odds scores (in bits) into p-values, assuming a zero-order background model. By default the program reports all motif occurrences with a p-value less than 1e-4. The threshold can be set using the --thresh option. The p-values for each motif occurence are converted to q-values following the method of Benjamini and Hochberg ("q-value" is defined as the minimal false discovery rate at which a given motif occurrence is deemed significant). The --qv-thresh option directs the program to use q-values rather than p-values for the threshold. If a motif has the strand feature set to +/- (rather than +), then fimo will search both strands for occurrences.

The parameter --max-stored-scores sets the maximum number of motif occurrences that will be retained in memory. It defaults to 100,000. If the number of matches found reaches the maximum value allowed, FIMO will discard 50% of the least significant matches, and new matches falling below the significance level of the retained matches will also be discarded.

FIMO can make use of position specific priors (PSP) to improve its identification of true motif occurrences. To take advantage of PSP in FIMO you use must provide two command line options. The --psp option is used to set the name of a MEME PSP file, and the --prior-dist option is used to set the name of a file containing the binned distribution of priors.


	USAGE: fimo [options] <motif file> <sequence file>

   	Options:
     --alpha <double> (default 1.0)
     --bgfile <background> (default from NR sequence database)
     --max-seq-length <int> (default=2.5e8)
     --max-stored-scores <int> (default=100000)
     --max-strand
     --motif <id> (default=all)
     --motif-pseudo <float> (default=0.1)
     --no-qvalue
     --norc
     --o <output dir> (default=fimo_out)
     --oc <output dir> (default=fimo_out)
     --parse-genomic-coord
     --psp <PSP filename> (default none)
     --prior-dist <PSP distribution filename> (default none)
     --qv-thresh
     --text
     --thresh <float> (default = 1e-4)
     --verbosity [1|2|3|4] (default 2)

   	Use '-' for <sequence file> to read the database from standard input.
   	Use '--bgfile motif-file' to read the background from the motif file.

--alpha <float> - The alpha parameter for calculating position specific priors. Alpha represents the fraction of all transcription factor binding sites that are binding sites for the TF of interest. Alpha must be between 0 and 1. The default value is 1.0.
--bgfile <bfile> - Read background frequencies from <bfile>. The file should be in MEME background file format. The default is to use frequencies embedded in the application from the non-redundant database. If the argument is the keyword motif-file, then the frequencies will be taken from the motif file.
--max-seq-length <max> - Set the maximum length allowed for input sequences. By default the maximum allowed length is 250000000.
--max-strand - If matches on both strands at a given position satisfy the output threshold, only report the match for the strand with the higher score. By default both matches are reported.
--max-stored-scores <max> - Set the maximum number of scores that will be stored. Precise calculation of q-values depends on having a complete list of scores. However, keeping a complete list of scores may exceed available memory. Once the number of stored scores reaches the maximum allowed, the least significant 50% of scores will be dropped, and approximate q-values will be calculated. By default the maximum number of stored matches is 100,000.
--motif <id> - Use only the motif identified by <id>. This option may be repeated.
--motif-pseudo <float> - A pseudocount to be added to each count in the motif matrix, after first multiplying by the corresponding background frequency (default=0.1).
--no-qvalue - Do not compute a q-value for each p-value. The q-value calculation is that of Benjamini and Hochberg (1995). By default, q-values are computed.
--norc - Do not score the reverse complement DNA strand. Both strands are scored by default.
--o <dir name> - Specifies the output directory. If the directory already exists, the contents will not be overwritten.
--oc <dir name> - Specifies the output directory. If the directory already exists, the contents will be overwritten.
--parse-genomic-coord When this option is specified each sequence header will be checked for UCSC style genomic coordinates. These are of the form >seq-name:starting pos.-ending pos.
where seq-name is the name of the sequence, starting pos. is the index of the first base, and ending pos. is the index of the final base. seq-name may not contain any white space. If genomic coordinates are found they will be used as the coordinates in the output. If no coordinate are found, the first position in the sequence will assumed to be 1.
--psp <file> - File containing position specific priors (PSP) in MEME PSP format.
--prior-dist <file> - File containing binned distribution of priors. This file can be generated from a MEME PSP format file. using the compute-prior-dist utility.
--qv-thresh - Directs the program to use q-values for the output threshold. The default is to use p-values.
--text Limits output to plain text sent to standard out. For FIMO, the text output is unsorted, and q-values are not reported. This mode allows the program to search an arbitrarily large database, because results are not stored in memory.
--thresh<float> - The output threshold for displaying search results. Only search results with a p-value less than the threshold will be output. The default threshold is a p-value of 1e-4. The threshold can be set to use q-values rather than p-values via the --qv-thresh option.
--verbosity 1|2|3|4 - Set the verbosity of status reports to standard error. The default level is 2.

	 */
	
	private String sequenceFasta;
	private String motifFile;
	private String outputDir;
	private String executable;
	private String additionalOptions;
	private String lsfJobID;
	private String description;
	private double alpha;
	private double qvalThreshold;
	
	/**
	 * Additional Fimo options
	 */
	public static String ADDITIONAL_OPTIONS = "-bgfile motif-file --norc --verbosity 4";
	
	
	/**
	 * @param executableFile Fimo executable
	 * @param sequenceFile Sequence fasta file
	 * @param motifXml Xml motif file
	 * @param outDir Output directory
	 * @param optionAlpha The fraction of sequences that are expected to harbor the motif
	 * @param optionQvalueThreshold Q value threshold
	 * @param addlOptions Additional Fimo options
	 * @param jobDescription Description
	 */
	public FimoJob(String executableFile, String sequenceFile, String motifXml, String outDir, double optionAlpha, double optionQvalueThreshold, String addlOptions, String jobDescription) {
		sequenceFasta = sequenceFile;
		executable = executableFile;
		motifFile = motifXml;
		outputDir = outDir;
		additionalOptions = addlOptions;
		description = jobDescription;
		alpha = optionAlpha;
		qvalThreshold = optionQvalueThreshold;
	}
	
	/**
	 * Submit the job to LSF
	 * @param queue Queue name
	 * @param memoryRequest Memory request in gigabytes
	 * @param scheduler 
	 * @return Job ID
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public Job submitJob(String queue, int memoryRequest, Scheduler scheduler) throws IOException, InterruptedException {
		lsfJobID = description + "_" + System.currentTimeMillis();
		String bsubOut = outputDir + "/" + lsfJobID + ".bsub";
		String command = executable + " ";
		command += "-oc " + outputDir + " ";
		command += "-alpha " + alpha + " ";
		command += "--qv-thresh --thresh " + qvalThreshold + " ";
		command += additionalOptions + " ";
		command += motifFile + " ";
		command += sequenceFasta;
		switch(scheduler) {
		case LSF:
			LSFJob job = new LSFJob(Runtime.getRuntime(), lsfJobID, command, bsubOut, queue, memoryRequest);
			job.submit();
			return job;
		default:
			throw new IllegalArgumentException("Scheduler " + scheduler.toString() + " not supported.");
		}
	}
	
	
}
