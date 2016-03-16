package util.programs.search;

import java.util.Collection;

import guttmanlab.core.util.MismatchGenerator;

/**
 * Search for imperfect matches
 * @author prussell
 *
 */
public class ImperfectKmerSearch extends PerfectKmerSearch {
	
	/**
	 * @param queries  Queries to look for in targets passed to search methods
	 * @param maxMismatches Max number of mismatches to report a match
	 */
	public ImperfectKmerSearch(Collection<String> queries, int maxMismatches) {
		super(MismatchGenerator.getRepresentatives(queries, maxMismatches));
	}
	
}
