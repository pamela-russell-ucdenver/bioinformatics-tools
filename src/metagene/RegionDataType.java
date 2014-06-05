package metagene;

import java.io.IOException;
import java.util.List;

import nextgen.core.annotation.Gene;

/**
 * @author prussell
 *
 */
public interface RegionDataType {
	
	/**
	 * Get the data for the region
	 * @param region Region
	 * @param reverseIfMinusOrientation Reverse data order if region is on minus strand
	 * @return Data values
	 * @throws IOException 
	 */
	public List<Double> getData(Gene region, boolean reverseIfMinusOrientation) throws IOException;
	
	/**
	 * Get the summary number for the region
	 * @param region Region
	 * @return Summary statistic
	 */
	public double getSummary(Gene region);
	
}
