package annotation;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;

import org.apache.log4j.Logger;

import nextgen.core.annotation.Annotation.Strand;
import nextgen.core.annotation.Gene;
import nextgen.core.utils.AnnotationUtils;

import broad.core.parser.CommandLineParser;
import broad.pda.annotation.BEDFileParser;

public class AnnotationExpander {
	
	private Map<String, Collection<Gene>> children;
	private Map<Gene, Gene> childToParent;
	private static Logger logger = Logger.getLogger(AnnotationExpander.class.getName());
	
	private AnnotationExpander(String childrenBed, String parentsBed) throws IOException {
		children = BEDFileParser.loadDataByChr(new File(childrenBed));
		Map<String, Collection<Gene>> parents = BEDFileParser.loadDataByChr(new File(parentsBed));
		childToParent = AnnotationUtils.mapChildToLargestParent(children, parents);
	}
	
	private Gene expandChild(Gene child, int upstreamDist, int downstreamDist) {
		if(!childToParent.containsKey(child)) {
			logger.warn("Skipping child " + child.getName() + " because does not have parent");
			return null;
		}
		try {
			Gene parent = childToParent.get(child);
			if(parent.getOrientation().equals(Strand.UNKNOWN)) {
				throw new IllegalArgumentException("Strand of parent must be known");
			}
			int startSign = parent.getOrientation().equals(Strand.POSITIVE) ? -1 : 1;
			int endSign = -1 * startSign;
			int newStart = parent.genomicToGenomicPositionWithOffset(child.getStart(), startSign * upstreamDist);
			int newEnd = parent.genomicToGenomicPositionWithOffset(child.getEnd(), endSign * downstreamDist);
			Gene rtrn = parent.trimAbsolute(newStart, newEnd);
			rtrn.setName(child.getName() + "_expanded_" + upstreamDist + "_" + downstreamDist);
			return rtrn;
		} catch (Exception e) {
			logger.warn("Caught exception on child " + child.getName()+ ". Returning null.");
			return null;
		}
	}
	
	private void expandAllChildrenAndWrite(int upstreamDist, int downstreamDist, String outputBed) throws IOException {
		logger.info("");
		logger.info("Writing expanded annotations to " + outputBed);
		FileWriter w = new FileWriter(outputBed);
		for(Gene child : childToParent.keySet()) {
			Gene expanded = expandChild(child, upstreamDist, downstreamDist);
			if(expanded != null) {
				w.write(expanded.toBED() + "\n");
			}
		}
		w.close();
		logger.info("Done writing file.");
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.setProgramDescription("Identify sub-annotations with their parent annotation and expand the ends of each child according to the gene model of the parent");
		p.addStringArg("-bc", "Bed file of child annotations", true);
		p.addStringArg("-bp", "Bed file of parent annotations", true);
		p.addStringArg("-o", "Output bed file of expanded children", true);
		p.addIntArg("-u", "Number of bases to expand in upstream direction", true);
		p.addIntArg("-d", "Number of bases to expand in downstream direction", true);
		p.parse(args);
		String childBed = p.getStringArg("-bc");
		String parentBed = p.getStringArg("-bp");
		String outputBed = p.getStringArg("-o");
		int upstream = p.getIntArg("-u");
		int downstream = p.getIntArg("-d");
		
		AnnotationExpander a = new AnnotationExpander(childBed, parentBed);
		a.expandAllChildrenAndWrite(upstream, downstream, outputBed);
		
		logger.info("");
		logger.info("All done.");

	}

}
