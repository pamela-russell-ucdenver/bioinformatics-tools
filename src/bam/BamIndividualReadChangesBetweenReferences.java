package bam;

import guttmanlab.core.util.CommandLineParser;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.HashSet;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.log4j.Logger;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMRecordQueryNameComparator;
import net.sf.samtools.util.CloseableIterator;

public class BamIndividualReadChangesBetweenReferences {
	
	
	/*
	 * Method:
	 * 2 bam files sorted by query name
	 * Iterate through both bam files concurrently
	 * For each query, sentify reference mappings in exactly one bam file
	 * For each query, construct all ordered pairs of references such that the read "moved" from first to second
	 * (reduntantly count every possible ordered pair)
	 * Have a special "reference" for no mapping
	 * Keep track of counts of ordered pairs of references
	 * Print table
	 */

	private static Logger logger = Logger.getLogger(BamIndividualReadChangesBetweenReferences.class.getName());
	private static final String NONE = "none";
	
	/**
	 * An ordered pair of strings
	 * @author prussell
	 *
	 */
	private static final class StringOrderedPair implements Comparable<StringOrderedPair> {
		
		private String s1;
		private String s2;
		
		/**
		 * @param s1 String 1
		 * @param s2 String 2
		 */
		public StringOrderedPair(String s1, String s2) {
			this.s1 = s1;
			this.s2 = s2;
		}
		
		public String toString() {
			return "(" + s1 + " -> " + s2 + ")";
		}

		public String getS1() {return s1;}
		public String getS2() {return s2;}
		
		public boolean equals(Object o) {
			if(!o.getClass().equals(StringOrderedPair.class)) return false;
			StringOrderedPair op = (StringOrderedPair)o;
			return s1.equals(op.getS1()) && s2.equals(op.getS2());
		}
		
		public int hashCode() {
			return toString().hashCode();
		}
		
		@Override
		public int compareTo(StringOrderedPair o) {
			int s1compare = s1.compareTo(o.getS1());
			if(s1compare != 0) return s1compare;
			return s2.compareTo(o.getS2());
		}
		
	}
	
	/**
	 * Count ordered pairs of targets not in both bam files for individual queries
	 * @param bam1 Bam file 1
	 * @param bam2 Bam file 2
	 * @return For each ordered pair, number of times a query has it (one member of pair in one bam file, missing other member in other bam file)
	 */
	private static Map<StringOrderedPair, Integer> targetShiftCounts(String bam1, String bam2) {
		logger.info("");
		logger.info("Getting target shift counts between " + bam1 + " and " + bam2 + "...");
		Map<StringOrderedPair, Integer> rtrn = new TreeMap<StringOrderedPair, Integer>();
		CloseableIterator<Collection<StringOrderedPair>> iter = new SymmetricDifferenceIterator(bam1, bam2);
		int numDone = 0;
		while(iter.hasNext()) {
			numDone++;
			if(numDone % 1000000 == 0) {
				logger.info("Finished " + numDone + " queries");
			}
			Collection<StringOrderedPair> pairs = iter.next();
			for(StringOrderedPair pair : pairs) {
				if(!rtrn.containsKey(pair)) {
					rtrn.put(pair, Integer.valueOf(0));
				}
				rtrn.put(pair, Integer.valueOf(rtrn.get(pair).intValue() + 1));
			}
		}
		iter.close();
		return rtrn;
	}
	
	
	private static void write(Map<StringOrderedPair, Integer> counts, String outFile, String rowHeader) throws IOException {
		FileWriter w = new FileWriter(outFile);
		for(StringOrderedPair pair : counts.keySet()) {
			w.write(rowHeader + "\t" + pair.toString() + "\t" + counts.get(pair) + "\n");
		}
		w.close();
	}
	
	
	/**
	 * A query and a set of targets
	 * Instances are mutable
	 * @author prussell
	 *
	 */
	private static final class QueryWithTargets {
		
		private String query;
		private Collection<String> targets;
		
		/**
		 * @param query The query
		 */
		public QueryWithTargets(String query) {
			this.query = query;
			this.targets = new HashSet<String>();
		}
		
		/**
		 * @param query The query
		 * @param target One target to initialize
		 */
		@SuppressWarnings("unused")
		public QueryWithTargets(String query, String target) {
			this.query = query;
			this.targets = new HashSet<String>();
			targets.add(target);
		}
		
		/**
		 * @param query The query
		 * @param targets Some targets to initialize
		 */
		@SuppressWarnings("unused")
		public QueryWithTargets(String query, Collection<String> targets) {
			this.query = query;
			this.targets = targets;
		}
		
		/**
		 * Add one target
		 * @param target Target to add
		 */
		public void addTarget(String target) {targets.add(target);}
		
		/**
		 * Get an object with a query name but no targets
		 * @param queryName Query name
		 * @return Object with query name and empty set of targets
		 */
		public static QueryWithTargets queryWithNoTargets(String queryName) {
			return new QueryWithTargets(queryName);
		}
		
		/**
		 * Add multiple targets
		 * @param targetsToAdd Targets to add
		 */
		@SuppressWarnings("unused")
		public void addTargets(Collection<String> targetsToAdd) {targets.addAll(targetsToAdd);}
		
		public String getQuery() {return query;}
		
		public Collection<String> getTargets() {return targets;}
		
	}
	
	/**
	 * Iterate over a bam file and get successive groups of records for the same query
	 * Bam file must be sorted by query name
	 * Skip unmapped records
	 * @author prussell
	 *
	 */
	private static final class QueryGroupIterator implements CloseableIterator<QueryWithTargets> {
		
		private SAMFileReader reader;
		private SAMRecordIterator samIter;
		private SAMRecord alreadyConsumedRecordThisCluster;
		private long numDone;
		
		public QueryGroupIterator(String bamFile) {
			numDone = 0;
			reader =  new SAMFileReader(new File(bamFile));
			samIter = reader.iterator();
			SAMFileHeader header = reader.getFileHeader();
			if(!header.getSortOrder().equals(SAMFileHeader.SortOrder.queryname)) {
				throw new IllegalArgumentException("Bam file must be sorted by query name");
			}
			//samIter.assertSorted(SAMFileHeader.SortOrder.queryname);
			if(!samIter.hasNext()) {
				throw new IllegalStateException("SAM iterator is empty");
			}
			alreadyConsumedRecordThisCluster = samIter.next();
		}
		
		@Override
		public boolean hasNext() {
			return samIter.hasNext();
		}

		@Override
		public QueryWithTargets next() {
			
			// Initialize with previously read record
			String query = alreadyConsumedRecordThisCluster.getReadName();
			QueryWithTargets rtrn = new QueryWithTargets(query);
			if(!alreadyConsumedRecordThisCluster.getReadUnmappedFlag()) {rtrn.addTarget(alreadyConsumedRecordThisCluster.getReferenceName());}
			
			while(samIter.hasNext()) {
				SAMRecord next = samIter.next();
				numDone++;
				if(numDone % 1000000 == 0) {
					logger.info("Finished " + numDone + " records");
				}
				if(next.getReadName().equals(query)) {
					if(next.getReadUnmappedFlag()) {
						continue;
					} else {rtrn.addTarget(next.getReferenceName());}
				} else {
					alreadyConsumedRecordThisCluster = next;
					break;
				}
			}
			return rtrn;
		}

		@Override
		public void close() {
			samIter.close();
			reader.close();
		}
		
		
	}
	
	/**
	 * Iterate through 2 query-sorted bam files and for each query, get all ordered pairs of targets in only one of the bams
	 * Ordered pair is a pair of one target only in bam file 1 and one target only in bam file 2
	 * If a query is only mapped in one bam file, ordered pairs will include the "none" string for the other bam file
	 * @author prussell
	 *
	 */
	private static final class SymmetricDifferenceIterator implements CloseableIterator<Collection<StringOrderedPair>> {

		private QueryGroupIterator groupIter1;
		private QueryGroupIterator groupIter2;
		private QueryWithTargets curr1;
		private QueryWithTargets curr2;
		
		/**
		 * @param bam1 Bam file 1
		 * @param bam2 Bam file 2
		 */
		public SymmetricDifferenceIterator(String bam1, String bam2) {
			groupIter1 = new QueryGroupIterator(bam1);
			groupIter2 = new QueryGroupIterator(bam2);
			if(!groupIter1.hasNext()) {
				throw new IllegalStateException("Iterator over " + bam1 + " is empty");
			}
			if(!groupIter2.hasNext()) {
				throw new IllegalStateException("Iterator over " + bam2 + " is empty");
			}
			curr1 = groupIter1.next();
			curr2 = groupIter2.next();
		}
		
		@Override
		public boolean hasNext() {
			return groupIter1.hasNext() || groupIter2.hasNext();
		}

		@Override
		public Collection<StringOrderedPair> next() {
			//logger.info("");
			if(curr1 == null) {
				//logger.info("Curr1 is null");
				String query = curr2.getQuery();
				//logger.info("Query from curr2 is " + query);
				Collection<StringOrderedPair> rtrn = getAllPairsOfTargetsNotInBoth(QueryWithTargets.queryWithNoTargets(query), curr2);
				//for(StringOrderedPair pair : rtrn) logger.info(pair.toString());
				curr2 = groupIter2.next();
				return rtrn;
			}
			if(curr2 == null) {
				//logger.info("Curr2 is null");
				String query = curr1.getQuery();
				//logger.info("Query from curr1 is " + query);
				Collection<StringOrderedPair> rtrn = getAllPairsOfTargetsNotInBoth(curr1, QueryWithTargets.queryWithNoTargets(query));
				//for(StringOrderedPair pair : rtrn) logger.info(pair.toString());
				curr1 = groupIter1.next();
				return rtrn;
			}
			int queryCompare = SAMRecordQueryNameComparator.compareReadNames(curr1.getQuery(), curr2.getQuery());
			if(queryCompare == 0) {
				//logger.info("Query 1 (" + curr1.getQuery() + ") is equal to query 2 (" + curr2.getQuery() + ")");
				Collection<StringOrderedPair> rtrn = getAllPairsOfTargetsNotInBoth(curr1, curr2);
				//for(StringOrderedPair pair : rtrn) logger.info(pair.toString());
				curr1 = groupIter1.next();
				curr2 = groupIter2.next();
				return rtrn;
			}
			if(queryCompare < 0) {
				//logger.info("Query 1 (" + curr1.getQuery() + ") is before query 2 (" + curr2.getQuery() + "). Using query 1 only.");
				String query = curr1.getQuery();
				Collection<StringOrderedPair> rtrn = getAllPairsOfTargetsNotInBoth(curr1, QueryWithTargets.queryWithNoTargets(query));
				//for(StringOrderedPair pair : rtrn) logger.info(pair.toString());
				curr1 = groupIter1.next();
				return rtrn;
			}
			//logger.info("Query 1 (" + curr1.getQuery() + ") is after query 2 (" + curr2.getQuery() + "). Using query 2 only.");
			String query = curr2.getQuery();
			Collection<StringOrderedPair> rtrn = getAllPairsOfTargetsNotInBoth(QueryWithTargets.queryWithNoTargets(query), curr2);
			//for(StringOrderedPair pair : rtrn) logger.info(pair.toString());
			curr2 = groupIter2.next();
			return rtrn;
		}

		@Override
		public void close() {
			groupIter1.close();
			groupIter2.close();
		}
		
		
		
	}
	
	
	
	/**
	 * @return True iff the two objects have the same query string
	 */
	private static boolean sameQuery(QueryWithTargets q1, QueryWithTargets q2) {
		return q1.getQuery().equals(q2.getQuery());
	}
	
	/**
	 * Get all ordered pairs of target sequence names in set 1 only and set 2 only
	 * Targets in both sets are ignored
	 * The pairs are ordered: target from set 1 first, set 2 second
	 * @param q1 Query/targets object 1
	 * @param q2 Query/targets object 2
	 * @return All possible ordered pairs of a target in set 1 only and a target in set 2 only,
	 * using the NONE string if one of these difference sets is empty
	 */
	private static Collection<StringOrderedPair> getAllPairsOfTargetsNotInBoth(QueryWithTargets q1, QueryWithTargets q2) {
		if(!sameQuery(q1, q2)) {
			throw new IllegalArgumentException("Both sets must have same query");
		}
		return getAllPairsNotInBoth(q1.getTargets(), q2.getTargets());
	}
	
	/**
	 * Get all ordered pairs of a string in set 1 only and a string in set 2 only
	 * Strings in both sets are ignored
	 * The pairs are ordered: element of set 1 first, element of set 2 second
	 * If lists are passed in, equal elements are collapsed
	 * @param strings1 Set 1 of strings
	 * @param strings2 Set 2 of strings
	 * @return All possible ordered pairs of an element of set 1 not in set 2 
	 * and an element of set 2 not in set 1, using the NONE string if one of these difference sets is empty
	 */
	private static Collection<StringOrderedPair> getAllPairsNotInBoth(Collection<String> strings1, Collection<String> strings2) {
		
		Collection<String> onlyIn1 = new HashSet<String>();
		Collection<String> onlyIn2 = new HashSet<String>();
		for(String s : strings1) {if(!strings2.contains(s)) {onlyIn1.add(s);}}
		for(String s : strings2) {if(!strings1.contains(s)) {onlyIn2.add(s);}}
		Collection<StringOrderedPair> rtrn = new TreeSet<StringOrderedPair>();
		if(onlyIn1.isEmpty() && onlyIn2.isEmpty()) {
			// The sets were the same. Return empty collection.
			return rtrn;
		}
		
		if(onlyIn1.isEmpty()) {onlyIn1.add(NONE);}
		if(onlyIn2.isEmpty()) {onlyIn2.add(NONE);}
		
		for(String s1 : onlyIn1) {
			for(String s2: onlyIn2) {
				rtrn.add(new StringOrderedPair(s1, s2));
			}
		}
		
		return rtrn;
		
	}
	
	
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-b1", "Bam 1", true);
		p.addStringArg("-b2", "Bam 2", true);
		p.addStringArg("-r", "Row header for every row of output table", true);
		p.addStringArg("-o", "Out table", true);
		p.parse(args);
		String bam1 = p.getStringArg("-b1");
		String bam2 = p.getStringArg("-b2");
		String out = p.getStringArg("-o");
		String rowHeader = p.getStringArg("-r");
		
		Map<StringOrderedPair, Integer> counts = targetShiftCounts(bam1, bam2);
		write(counts, out, rowHeader);
		
		logger.info("");
		logger.info("All done.");
		
	}
	
	
}
