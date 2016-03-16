package util.programs.counts;

import guttmanlab.core.util.CommandLineParser;
import guttmanlab.core.util.StringParser;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.apache.log4j.Logger;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFormatException;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class BamFileCountsByReferenceCategory {
	
	private Map<String, Category> categoriesByName;
	private Map<String, String> categoryNamesBySeqName;
	private static Logger logger = Logger.getLogger(BamFileCountsByReferenceCategory.class.getName());
	private Category categoryOther;
	private Category categoryUnmapped;

	private class Category {
		
		private String name;
		private int count;
		
		public Category(String categoryName) {
			count = 0;
			name = categoryName;
		}
		
		public void incrementCount() {
			count++;
		}
		
		public int getCount() {
			return count;
		}
		
		public String getName() {
			return name;
		}
		
	}
	
	private BamFileCountsByReferenceCategory(String bamFile, String categoryTable) throws IOException {
		categoryOther = new Category("other");
		categoryUnmapped = new Category("unmapped");
		loadCategoryTable(categoryTable);
		makeCounts(bamFile);
	}
	
	private void loadCategoryTable(String tableFile) throws IOException {
		logger.info("");
		categoriesByName = new HashMap<String, Category>();
		categoryNamesBySeqName = new HashMap<String, String>();
		BufferedReader r = new BufferedReader(new FileReader(tableFile));
		StringParser s = new StringParser();
		while(r.ready()) {
			s.parse(r.readLine());
			if(s.getFieldCount() == 0) continue;
			if(s.getFieldCount() != 2) {
				r.close();
				throw new IllegalArgumentException("Illegal number of fields: " + s.getFieldCount() + ". Table format: sequence_name   category_name");
			}
			String seqName = s.asString(0);
			String catName = s.asString(1);
			if(categoryNamesBySeqName.containsKey(seqName)) {
				r.close();
				throw new IllegalStateException("Sequence name " + seqName + " is duplicated");
			}
			if(!categoriesByName.containsKey(catName)) {
				categoriesByName.put(catName, new Category(catName));
				logger.info("Added category " + categoriesByName.get(catName).getName());
			}
			categoryNamesBySeqName.put(seqName, catName);
		}
		r.close();
	}
		
	private void makeCounts(String bamFile) {
		logger.info("");
		logger.info("Making counts...");
		SAMFileReader reader = new SAMFileReader(new File(bamFile));
		SAMRecordIterator iter = reader.iterator();
		int numDone = 0;
		
		while(iter.hasNext()) {
			try {
				SAMRecord record = iter.next();
				numDone++;
				if(numDone % 1000000 == 0) {
					logger.info("Finished " + numDone + " records.");
				}
				if(!record.getReadUnmappedFlag()) {
					String refName = record.getReferenceName();
					if(!categoryNamesBySeqName.containsKey(refName)) {
						categoryOther.incrementCount();
						continue;
					}
					String categoryName = categoryNamesBySeqName.get(refName);
					categoriesByName.get(categoryName).incrementCount();
				} else {
					categoryUnmapped.incrementCount();
				}
			} catch(SAMFormatException e) {
				logger.info("Skipping record: " + e.getMessage());
			}
		}
		
		reader.close();
		logger.info("Done making counts.");
	}
	
	private void printCounts(String outFile) throws IOException {
		FileWriter w = new FileWriter(outFile);
		for(Category category : categoriesByName.values()) {
			w.write(category.getName() + "\t" + category.getCount() + "\n");
		}
		w.write(categoryOther.getName() + "\t" + categoryOther.getCount() + "\n");
		w.write(categoryUnmapped.getName() + "\t" + categoryUnmapped.getCount() + "\n");
		w.close();
	}
	
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-t", "Table of <reference_name>  <category_name>", true);
		p.addStringArg("-b", "Bam file", true);
		p.addStringArg("-o", "Output table of counts", true);
		p.parse(args);
		String categoryTable = p.getStringArg("-t");
		String bamFile = p.getStringArg("-b");
		String outFile = p.getStringArg("-o");
		
		BamFileCountsByReferenceCategory b = new BamFileCountsByReferenceCategory(bamFile, categoryTable);
		b.printCounts(outFile);

	}

}
