package util.programs.bam;

import guttmanlab.core.util.CommandLineParser;
import htsjdk.samtools.fork.util.iterators.MappingCountPerQuery;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Average number of mappings per read
 * @author prussell
 *
 */
public class MappingsPerQuery {
	
	private static double meanCountsPerQuery(String bamFile) {
		List<Integer> countsList = new ArrayList<Integer>();
		MappingCountPerQuery iter = new MappingCountPerQuery(bamFile);
		while(iter.hasNext()) countsList.add(iter.next());
		iter.close();
		return countsList.stream().collect(Collectors.summarizingInt(x -> x)).getAverage();
	}
	
	public static void main(String[] args) {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-b", "Bam file", true);
		p.parse(args);
		String bam = p.getStringArg("-b");
		double mean = meanCountsPerQuery(bam);
		System.out.println(bam + "\t" + mean);
		
	}
	
	
}
