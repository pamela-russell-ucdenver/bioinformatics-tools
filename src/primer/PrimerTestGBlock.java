package primer;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;

import org.apache.log4j.Logger;

import broad.core.parser.CommandLineParser;
import broad.core.parser.StringParser;
import broad.core.primer3.PrimerPair;
import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;

public class PrimerTestGBlock {
	
	private Collection<PrimerPair> primers;
	private int primerLength;
	private static Logger logger = Logger.getLogger(PrimerTestGBlock.class.getName());
	
	private PrimerTestGBlock(String primerFile) throws IOException {
		primers = new ArrayList<PrimerPair>();
		FileReader r = new FileReader(primerFile);
		BufferedReader b = new BufferedReader(r);
		StringParser s = new StringParser();
		while(b.ready()) {
			s.parse(b.readLine());
			primers.add(new PrimerPair(s.getStringArray()));
		}
		b.close();
		primerLength = primers.iterator().next().getLeftPrimer().length();
		for(PrimerPair primer : primers) {
			if(primer.getLeftPrimer().length() != primerLength || primer.getRightPrimer().length() != primerLength) {
				throw new IllegalArgumentException("All primers must have same length");
			}
		}
	}
	
	private List<Sequence> createGBlocks(int numPrimersPerSegment, int numSegmentsPerBlock) {
		
		int numPrimersPerBlock = numPrimersPerSegment * numSegmentsPerBlock;
		int blockLength = primerLength * 2 * numPrimersPerBlock;
		int numBlocks = (primers.size() / numPrimersPerBlock);
		if(primers.size() % numPrimersPerBlock != 0) {
			numBlocks++;
		}
		logger.info("Creating " + numBlocks + " blocks of length " + blockLength + " with " + numPrimersPerSegment + " primers per segment and " + numPrimersPerBlock + " primers per block.");
		if(primers.size() % numPrimersPerBlock != 0) {
			throw new IllegalArgumentException("Total number of primers must be divisible by " + numPrimersPerBlock);
		}
		
		List<Sequence> rtrn = new ArrayList<Sequence>();
		
		Iterator<PrimerPair> primerIter = primers.iterator();
		
		for(int currBlock = 0; currBlock < numBlocks; currBlock++) {
			String[] block = new String[numPrimersPerBlock * 2];
			for(int currSegment = 0; currSegment < numSegmentsPerBlock; currSegment++) {
				for(int currPrimer = 0; currPrimer < numPrimersPerSegment; currPrimer++) {
					PrimerPair primer = primerIter.next();
					int left = currSegment * 2 * numPrimersPerSegment + currPrimer;
					int right = numPrimersPerSegment + currSegment * 2 * numPrimersPerSegment + currPrimer;
					block[left] = primer.getLeftPrimer();
					block[right] = Sequence.reverseSequence(primer.getRightPrimer());
					//logger.info("Block " + currBlock + " segment " + currSegment + " primer " + currPrimer + " left " + left + " right " + right);
				}
			}
			//logger.info("");
			//logger.info("Block " + currBlock);
			for(int i = 0; i < block.length; i++) {
				//logger.info(i + "\t" + block[i]);
			}
			String blockString = "";
			for(int i = 0; i < block.length; i++) {
				blockString += block[i];
			}
			Sequence blockSeq = new Sequence("block_" + currBlock);
			blockSeq.setSequenceBases(blockString);
			rtrn.add(blockSeq);
		}
		
		return rtrn;
		
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-p", "Primer file", true);
		p.addIntArg("-ps", "Primers per segment", true);
		p.addIntArg("-sb", "Segments per block", true);
		p.addStringArg("-o", "Output file", true);
		p.parse(args);
		String primerFile = p.getStringArg("-p");
		int primersPerSegment = p.getIntArg("-ps");
		int segmentsPerBlock = p.getIntArg("-sb");
		String outFile = p.getStringArg("-o");
		
		PrimerTestGBlock b = new PrimerTestGBlock(primerFile);
		List<Sequence> gblocks = b.createGBlocks(primersPerSegment, segmentsPerBlock);
		
		FastaSequenceIO fsio = new FastaSequenceIO(outFile);
		fsio.write(gblocks);
		

	}

}
