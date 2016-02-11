package bam;

import guttmanlab.core.util.CommandLineParser;
import guttmanlab.core.util.StringParser;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.apache.log4j.Logger;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class MultimapperCountsToEquivalentTargets {
	
	private Map<String, String> targetToEquivClass;
	private SAMFileReader samReader;
	private SAMRecordIterator samRecordIter;
	private Map<QueryMappingStatus, Integer> mappingStatusCounts;
	private SAMRecord currentRecord;
	private static Logger logger = Logger.getLogger(MultimapperCountsToEquivalentTargets.class.getName());
	private static boolean plusStrandOnly = false;
	private int numRecordsDone;
	
	private enum QueryMappingStatus {
		UNMAPPED,
		SINGLE_EQUIV_CLASS,
		MULTIPLE_EQUIV_CLASSES;
	}
	
	private MultimapperCountsToEquivalentTargets(String bamFile, String targetEquivClassTable) throws IOException {
		
		numRecordsDone = 0;
		
		// Bam file must be query sorted
		samReader = new SAMFileReader(new File(bamFile));
		SAMFileHeader header = samReader.getFileHeader();
		if(header.getSortOrder() != SAMFileHeader.SortOrder.queryname) {
			throw new IllegalArgumentException("Bam file must be sorted by query name");
		}
		samRecordIter = samReader.iterator();
		if(samRecordIter.hasNext()) {
			currentRecord = samRecordIter.next();
		}
		
		BufferedReader b = new BufferedReader(new FileReader(targetEquivClassTable));
		targetToEquivClass = new HashMap<String, String>();
		StringParser s = new StringParser();
		while(b.ready()) {
			s.parse(b.readLine());
			if(s.getFieldCount() != 2) {
				b.close();
				throw new IllegalArgumentException("Line format: target_name   equivalence_class");
			}
			targetToEquivClass.put(s.asString(0), s.asString(1));
		}
		b.close();
	}
	
	
	private List<SAMRecord> nextQuery() {
		List<SAMRecord> rtrn = new ArrayList<SAMRecord>();
		while(samRecordIter.hasNext()) {
			SAMRecord next = samRecordIter.next();
			numRecordsDone++;
			if(numRecordsDone % 1000000 == 0) {
				logger.info("Finished " + numRecordsDone + " records.");
			}
			if(plusStrandOnly) {
				if(next.getReadNegativeStrandFlag()) {
					continue;
				}
			}
			if(next.getReadName().equals(currentRecord.getReadName())) {
				rtrn.add(next);
			}
			else {
				rtrn.add(currentRecord);
				currentRecord = next;
				return rtrn;
			}
		}
		return rtrn;
	}
	
	private QueryMappingStatus mappingStatus(Collection<SAMRecord> mappingsForQuery) {
		String query = mappingsForQuery.iterator().next().getReadName();
		Collection<String> mappedEquivClasses = new HashSet<String>();
		for(SAMRecord record : mappingsForQuery) {
			if(record.getReadUnmappedFlag()) {
				continue;
			}
			String readName = record.getReadName();
			if(!readName.equals(query)) {
				throw new IllegalArgumentException("Records don't have same query name: " + query + " " + record.getReadName());
			}
			mappedEquivClasses.add(targetToEquivClass.get(record.getReferenceName()));
		}
		switch(mappedEquivClasses.size()) {
		case 0:
			return QueryMappingStatus.UNMAPPED;
		case 1:
			return QueryMappingStatus.SINGLE_EQUIV_CLASS;
		default:
			return QueryMappingStatus.MULTIPLE_EQUIV_CLASSES;
		}
	}
	
	private void incrementStatusCount(QueryMappingStatus status) {
		mappingStatusCounts.put(status, Integer.valueOf(mappingStatusCounts.get(status).intValue() + 1));
	}
	
	private void countMappingStatuses() {
		mappingStatusCounts = new HashMap<QueryMappingStatus, Integer>();
		for(int i = 0; i < QueryMappingStatus.values().length; i++) {
			mappingStatusCounts.put(QueryMappingStatus.values()[i], Integer.valueOf(0));
		}
		while(true) {
			List<SAMRecord> nextQueryMappings = nextQuery();
			if(nextQueryMappings.isEmpty()) {
				break;
			}
			incrementStatusCount(mappingStatus(nextQueryMappings));
		}
		String message = "SINGLE_EQUIV_CLASS:" + mappingStatusCounts.get(QueryMappingStatus.SINGLE_EQUIV_CLASS).toString() + "\t";
		message += "MULTIPLE_EQUIV_CLASSES:" + mappingStatusCounts.get(QueryMappingStatus.MULTIPLE_EQUIV_CLASSES).toString() + "\t";
		message += "UNMAPPED:" + mappingStatusCounts.get(QueryMappingStatus.UNMAPPED).toString() + "\t";
		logger.info(message);
	}
	
	
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-b", "Bam file", true);
		p.addStringArg("-e", "Equivalence class table for target sequences. Format: sequence   equivalence_class", true);
		p.addBooleanArg("-p", "Include plus strand mappings only", false, false);
		p.parse(args);
		String bamFile = p.getStringArg("-b");
		String equivClassFile = p.getStringArg("-e");
		plusStrandOnly = p.getBooleanArg("-p");
		
		MultimapperCountsToEquivalentTargets m = new MultimapperCountsToEquivalentTargets(bamFile, equivClassFile);
		m.countMappingStatuses();
		
		logger.info("");
		logger.info("All done.");
		
	}
	
	
}
