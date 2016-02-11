package util;

import java.util.Collection;

/**
 * Search target string for one of a collection of queries stored in a hash
 * @author prussell
 *
 */
public class PerfectKmerSearch {
	
	private Collection<String> queries;
	private int queryLength;
	
	/**
	 * @param queryStrings Queries to look for in targets passed to search methods
	 */
	public PerfectKmerSearch(Collection<String> queryStrings) {
		if(queryStrings.isEmpty()) {
			throw new IllegalArgumentException("Must provide at least one query string");
		}
		queryLength = queryStrings.iterator().next().length();
		for(String s : queryStrings) {
			if(s.length() != queryLength) {
				throw new IllegalArgumentException("All query strings must have same length");
			}
		}
		queries = queryStrings;
	}
	
	/**
	 * Check if the target contains some query
	 * @param target Target string
	 * @return True iff target contains a perfect match for at least one of the queries
	 */
	public boolean containsQuery(String target) {
		for(int i = 0; i <= target.length() - queryLength; i++) {
			if(queries.contains(target.substring(i, i + queryLength))) {
				return true;
			}
		}
		return false;
	}
	
}
