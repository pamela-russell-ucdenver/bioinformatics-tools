package bed;

import general.CommandLineParser;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.log4j.Logger;

import nextgen.core.annotation.Gene;

import broad.pda.annotation.BEDFileParser;

public class BedCollapseSameNamesAndEnds {
	
	private static int dist(int a, int b) {
		return Math.max(a-b, b-a);
	}
	
	private static boolean areClose(Gene gene1, Gene gene2, int max) {
		return dist(gene1.getStart(), gene2.getStart()) <= max && dist(gene1.getEnd(), gene2.getEnd()) <= max;
	}
	
	private static Logger logger = Logger.getLogger(BedCollapseSameNamesAndEnds.class.getName());
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-i", "Input bed", true);
		p.addStringArg("-o", "Output bed", true);
		p.addIntArg("-d", "Max distance to collapse endpoints", true);
		p.parse(args);
		String input = p.getStringArg("-i");
		String output = p.getStringArg("-o");
		int dist = p.getIntArg("-d");
		
		Map<String, Collection<Gene>> genesByChr = BEDFileParser.loadDataByChr(new File(input));
		Map<String, TreeSet<Gene>> genesByName = new TreeMap<String, TreeSet<Gene>>();
		for(String chr : genesByChr.keySet()) {
			for(Gene gene : genesByChr.get(chr)) {
				String name = gene.getName();
				if(!genesByName.containsKey(name)) {
					genesByName.put(name, new TreeSet<Gene>());
				}
				genesByName.get(name).add(gene);
			}
		}
		
		for(String name : genesByName.keySet()) {
			Iterator<Gene> iter = genesByName.get(name).iterator();
			Gene lastKept = iter.next();
			while(iter.hasNext()) {
				Gene next = iter.next();
				if(areClose(lastKept, next, dist)) {
					logger.info("Removing " + next.getName() + ":" + next.toUCSC() + " because it is too close to " + lastKept.getName() + ":" + lastKept.toUCSC());
					iter.remove();
					continue;
				}
				lastKept = next;
			}
		}
		
		FileWriter w = new FileWriter(output);
		for(String name : genesByName.keySet()) {
			int i = 0;
			for(Gene gene : genesByName.get(name)) {
				gene.setName(name + "_endpoints_group_" + i);
				i++;
				w.write(gene.toBED() + "\n");
			}
		}
		w.close();
		
		logger.info("");
		logger.info("All done.");
		
	}

}
