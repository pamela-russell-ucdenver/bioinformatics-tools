package metagene;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.NavigableMap;
import java.util.TreeMap;


/**
 * @author prussell
 *
 */
public final class Util {

	
	/**
	 * Get space delimited string of all elements in list
	 * @param list The list
	 * @return String representation of list
	 */
	public static String listToString(List<Double> list) {
		String rtrn = "";
		for(Double d : list) {
			rtrn += d.toString() + "\t";
		}
		return rtrn;
	}
	
	/**
	 * Get table representing the list
	 * First column is list position; second column is value
	 * @param list List of values
	 * @return Table as string with newline characters
	 */
	public static String listToTableByPos(List<Double> list) {
		String rtrn = "";
		for(int i=0; i<list.size(); i++) {
			rtrn += i + "\t" + list.get(i).toString() + "\n";
		}
		return rtrn;
	}
	
	/**
	 * Sum two lists of equal size
	 * @param list1 List 1
	 * @param list2 List 2
	 * @return List of the sum by position
	 */
	public static List<Double> sum(List<Double> list1, List<Double> list2) {
		if(list1.size() != list2.size()) {
			throw new IllegalArgumentException("Lists must have same size.");
		}
		List<Double> rtrn = new ArrayList<Double>();
		for(int i = 0; i < list1.size(); i++) {
			rtrn.add(Double.valueOf(list1.get(i).doubleValue() + list2.get(i).doubleValue()));
		}
		return rtrn;
	}
	
	/**
	 * Sum a collection of lists of equal size
	 * @param lists The lists
	 * @return List of the sum by position
	 */
	public static List<Double> sum(Collection<List<Double>> lists) {
		int size = lists.iterator().next().size();
		for(List<Double> list : lists) {
			if(list.size() != size) {
				throw new IllegalArgumentException("All lists must have same size.");
			}
		}
		List<Double> rtrn = new ArrayList<Double>();
		for(int i=0; i<size; i++) {
			double sum = 0;
			for(List<Double> list : lists) {
				sum += list.get(i).doubleValue();
			}
			rtrn.add(Double.valueOf(sum));
		}
		return rtrn;
	}
	
	/**
	 * Get the sum of numbers in a list
	 * @param list The list
	 * @return The sum
	 */
	public static double sum(List<Double> list) {
		double rtrn = 0;
		for(Double d : list) {
			rtrn += d.doubleValue();
		}
		return rtrn;
	}
	
	/**
	 * Scale a list to control the sum of elements
	 * @param list The list
	 * @param newSum The desired sum
	 * @return Normalized list
	 */
	public static List<Double> normalizeSum(List<Double> list, double newSum) {
		double sum = sum(list);
		double normalizationFactor = newSum / sum;
		List<Double> rtrn = new ArrayList<Double>();
		for(Double d : list) {
			rtrn.add(Double.valueOf(d.doubleValue() * normalizationFactor));
		}
		return rtrn;
	}
	
	/**
	 * Contract or expand a list to a specified size by averaging over multiple entries
	 * @param list The list to scale
	 * @param resultSize The desired size
	 * @return New list of the desired size that represents a smoothing of the original list
	 */
	public static List<Double> expandOrContractList(List<Double> list, int resultSize) {
		
		if(list.size() == resultSize) {
			return list;
		}
		
		
		// Give each list element a normalized position relative to the result size
		TreeMap<Double, Double> listWithNormalizedPositions = new TreeMap<Double,Double>();
		for(int i=0; i < list.size(); i++) {
			double normalizedPos = resultSize * ((double)i / (double)list.size());
			listWithNormalizedPositions.put(Double.valueOf(normalizedPos), list.get(i));
		}
		
		List<Double> rtrn = new ArrayList<Double>();
		for(int i=0; i<resultSize; i++) {
			
			// Try to get all elements whose normalized position is between i and i+1
			NavigableMap<Double, Double> subtree = listWithNormalizedPositions.subMap(Double.valueOf(i), true, Double.valueOf(i+1), true);
			
			// If there is nothing in the interval, average the previous and next elements
			if(subtree.size() == 0) {
				
				Double avg = null;

				// Get the previous and next keys
				Double previousKey = listWithNormalizedPositions.floorKey(Double.valueOf(i));
				Double nextKey = listWithNormalizedPositions.ceilingKey(Double.valueOf(i));
				if(previousKey == null && nextKey == null) {
					throw new IllegalStateException("No previous or next key for normalized position " + i);
				}
				if(previousKey == null) {
					// If there is no previous key use the next key
					avg = listWithNormalizedPositions.get(nextKey);
				} else if(nextKey == null) {
					// If there is no next key use the previous key
					avg = listWithNormalizedPositions.get(previousKey);
				} else {
					// If there are both keys use the average of previous and next
					double previous = listWithNormalizedPositions.get(previousKey).doubleValue();
					double next = listWithNormalizedPositions.get(nextKey).doubleValue();
					avg = Double.valueOf((previous + next) / 2);
				}
				rtrn.add(avg);
				
			} else {
				// If there are elements with normalized position between i and i+1, average all of them
				double total = 0;
				for(Double value : subtree.values()) {
					total += value.doubleValue();
				}
				double avg = total / subtree.size();
				rtrn.add(Double.valueOf(avg));
			}
		}
		return rtrn;
	}
	

}
