package util.programs.text;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;

import util.FileUtils;
import guttmanlab.core.util.CommandLineParser;

/**
 * Combine a set of strings into equivalence classes based on which strings contain common substrings ("markers")
 * @author prussell
 *
 */
public class CombineStringsByEquivalence {
	
	private static Logger logger = Logger.getLogger(CombineStringsByEquivalence.class.getName());
	
	private CombineStringsByEquivalence() {}
	
	/**
	 * A collection of strings that are considered to be equivalent
	 * And the set of markers (substrings contained in the strings) that make them equivalent
	 * @author prussell
	 *
	 */
	private class EquivalentStrings {
		
		private Set<String> strings;
		private Set<String> markers;
		
		public EquivalentStrings() {
			strings = new HashSet<String>();
			markers = new HashSet<String>();
		}
		
		@SuppressWarnings("unused")
		public void addString(String s) {
			strings.add(s);
		}
		
		public void addAllStrings(Collection<String> s) {
			strings.addAll(s);
		}
		
		public void addMarker(String marker) {
			markers.add(marker);
		}
		
		@SuppressWarnings("unused")
		public void clear() {
			strings.clear();
			markers.clear();
		}
		
		@SuppressWarnings("unused")
		public void combine(EquivalentStrings other) {
			strings.addAll(other.getStrings());
			markers.addAll(other.getMarkers());
		}
		
		public Collection<String> getStrings() {
			return strings;
		}
		
		public Collection<String> getMarkers() {
			return markers;
		}
		
		public boolean containsMarker(String marker) {
			return markers.contains(marker);
		}
		
		public boolean containsString(String s) {
			return strings.contains(s);
		}
		
		/**
		 * Returns a combined list of the strings separated by a delimiter
		 */
		public String toString() {
			if(strings.isEmpty()) {
				return "";
			}
			Iterator<String> iter = strings.iterator();
			String rtrn = "";
			if(iter.hasNext()) {
				rtrn += iter.next();
			}
			while(iter.hasNext()) {
				rtrn += "%" + iter.next();
			}
			return rtrn;
		}
		
		public boolean equals(Object o) {
			if(!o.getClass().equals(getClass())) {
				return false;
			}
			EquivalentStrings es = (EquivalentStrings)o;
			return es.toString().equals(toString());
		}
		
		public int hashCode() {
			return toString().hashCode();
		}
		
	}
	
	/**
	 * Map each marker to the collection of strings that contain it as a substring
	 * Different value maps can contain the same string if a string includes multiple markers
	 * @param markers The markers
	 * @param strings The strings
	 * @param keepStringsWithNoMarkers Include strings that don't contain any markers as their own representatives
	 * @return Map of marker to strings that contain it
	 */
	private static Map<String, Collection<String>> markerToAllStrings(Collection<String> markers, Collection<String> strings, boolean keepStringsWithNoMarkers) {
		Map<String, Collection<String>> rtrn = new HashMap<String, Collection<String>>();
		for(String marker : markers) {
			rtrn.put(marker, new HashSet<String>());
		}
		for(String string : strings) {
			boolean hasMarker = false;
			for(String marker : rtrn.keySet()) {
				if(string.contains(marker)) {
					hasMarker = true;
					rtrn.get(marker).add(string);
				}
			}
			if(keepStringsWithNoMarkers && !hasMarker) {
				Collection<String> thisString = new HashSet<String>();
				thisString.add(string);
				rtrn.put(string, thisString);
			}
		}
		return rtrn;
	}
	
	/**
	 * Make sure each marker and string appears in only one EquivalentStrings object
	 * @param equivalentStrings Set of EquivalentStrings objects to check for duplicates
	 * @param markerToAllStrings Map of marker to all strings that contain the marker
	 */
	private static void validate(Set<EquivalentStrings> equivalentStrings, Map<String, Collection<String>> markerToAllStrings) {
		Set<String> markers = new HashSet<String>();
		Set<String> strings = new HashSet<String>();
		for(EquivalentStrings es : equivalentStrings) {
			for(String marker : es.getMarkers()) {
				if(markers.contains(marker)) {
					throw new IllegalStateException("Two sets contain marker " + marker);
				}
				markers.add(marker);
			}
			for(String string : es.getStrings()) {
				if(strings.contains(string)) {
					throw new IllegalStateException("Two sets contain string " + string);
				}
				strings.add(string);
			}
		}
		for(String marker : markerToAllStrings.keySet()) {
			for(String string : markerToAllStrings.get(marker)) {
				boolean containsString = false;
				for(EquivalentStrings es : equivalentStrings) {
					if(es.containsString(string)) {
						containsString = true;
					}
				}
				if(!containsString) {
					throw new IllegalStateException("No set contains string " + string);
				}
			}
		}
	}
	
	/**
	 * Make a set of EquivalentStrings objects such that each marker and string appears in only one
	 * @param markerToAllStrings Map of marker to all strings that contain the marker
	 * @return Set of EquivalentStrings objects
	 */
	private static Set<EquivalentStrings> equivalentStrings(Map<String, Collection<String>> markerToAllStrings) {
		
		Set<EquivalentStrings> rtrn = new HashSet<EquivalentStrings>();
		
		for(String marker : markerToAllStrings.keySet()) {
			boolean found = false;
			for(EquivalentStrings es : rtrn) {
				if(es.containsMarker(marker)) {
					es.addAllStrings(es.getStrings());
					found = true;
					break;
				}
			}
			if(found) continue;
			for(String string : markerToAllStrings.get(marker)) {
				for(EquivalentStrings es : rtrn) {
					if(es.containsString(string)) {
						es.addMarker(marker);
						es.addAllStrings(markerToAllStrings.get(marker));
						found = true;
						break;
					}
				}
				if(found) break;
			}
			if(!found) {
				EquivalentStrings es = new CombineStringsByEquivalence().new EquivalentStrings();
				es.addMarker(marker);
				es.addAllStrings(markerToAllStrings.get(marker));
				rtrn.add(es);
			}
		}
		
		validate(rtrn, markerToAllStrings);
		return rtrn;
		
	}
	
	/**
	 * Create map of string to a string representation of its equivalence class (list of all equivalent strings)
	 * @param markers Set of markers
	 * @param strings Set of strings
	 * @param keepStringsWithNoMarkers Include strings that don't contain any markers as their own representatives
	 * @return Map of string to list of equivalent strings based on containing common markers
	 */
	private static Map<String, String> stringToEquivalenceClass(Collection<String> markers, Collection<String> strings, boolean keepStringsWithNoMarkers) {
		Map<String, Collection<String>> markerToAllStrings = markerToAllStrings(markers, strings, keepStringsWithNoMarkers);
		Set<EquivalentStrings> equivalentStrings = equivalentStrings(markerToAllStrings);
		Map<String, String> rtrn = new HashMap<String, String>();
		for(EquivalentStrings es : equivalentStrings) {
			String combined = es.toString();
			for(String string : es.getStrings()) {
				rtrn.put(string, combined);
			}
		}
		return rtrn;
	}
	
	/**
	 * Write table of map of string to equivalence class
	 * @param markerListFile File containing list of markers
	 * @param stringListFile File containing list of strings
	 * @param keepStringsWithNoMarkers Include strings that don't contain any markers as their own representatives
	 * @param outFile Output file
	 * @throws IOException
	 */
	private static void writeTableStringToEquivalenceClass(String markerListFile, String stringListFile, String outFile, boolean keepStringsWithNoMarkers) throws IOException {
		Collection<String> markers = FileUtils.fileLinesAsList(markerListFile);
		Collection<String> strings = FileUtils.fileLinesAsList(stringListFile);
		Map<String, String> stringToEquivClass = stringToEquivalenceClass(markers, strings, keepStringsWithNoMarkers);
		FileWriter w = new FileWriter(outFile);
		for(String s : stringToEquivClass.keySet()) {
			w.write(s + "\t" + stringToEquivClass.get(s) + "\n");
		}
		w.close();
	}
	
	
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-m", "List of markers (substrings of strings) whose presence makes the strings equivalent", true);
		p.addStringArg("-s", "List of strings that are equivalent if they contain the same markers", true);
		p.addStringArg("-o", "Output table of string and equivalence class of strings", true);
		p.addBooleanArg("-n", "Include strings with no markers as their own representatives", false, true);
		p.parse(args);
		String markers = p.getStringArg("-m");
		String strings = p.getStringArg("-s");
		String out = p.getStringArg("-o");
		boolean includeNoMarker = p.getBooleanArg("-n");
		
		writeTableStringToEquivalenceClass(markers, strings, out, includeNoMarker);
		
		logger.info("");
		logger.info("All done.");
		
	}
	
	
	
	
}
