package tests;

import org.apache.log4j.Logger;

import broad.core.parser.CommandLineParser;

public class Sleep {
	
	private static long MILLIS_PER_MIN = 60000;
	private static Logger logger = Logger.getLogger(Sleep.class.getName());

	/**
	 * @param args
	 * @throws InterruptedException 
	 */
	public static void main(String[] args) throws InterruptedException {
		
		CommandLineParser p = new CommandLineParser();
		p.addIntArg("-m", "Number of minutes to run", true);
		p.parse(args);
		int mins = p.getIntArg("-m");
		
		if(mins < 1) {
			throw new IllegalArgumentException("Must provide minutes >= 1");
		}
		
		logger.info("Running for " + mins + " minutes...");
		
		for(int i = 1; i <= mins; i++) {
			Thread.sleep(MILLIS_PER_MIN);
			logger.info(i + " minutes done");
		}
		
		logger.info("All done.");
		
	}

}
