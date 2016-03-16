package test.berkeleydb;

import java.util.Map;

/**
 * An item
 * @author prussell
 *
 */
public interface Item {
	
	/**
	 * @return Item name
	 */
	public String getName();
	
	/**
	 * @return Item color
	 */
	public String getColor();
	
	/**
	 * @param store Store
	 * @return Price at the store
	 */
	public float getPrice(Store store);
	
	/**
	 * Add a vendor
	 * @param store Store that sells this
	 * @param price Price at the store
	 */
	public void addVendor(Store store, float price);
	
	/**
	 * Get all vendors that sell this
	 * @return Vendors and prices
	 */
	public Map<Store, Float> getVendors();
	
}
