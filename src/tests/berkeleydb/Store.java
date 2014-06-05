package tests.berkeleydb;


/*
 * http://docs.oracle.com/cd/E17277_02/html/GettingStartedGuide/persistobject.html
 */
/**
 * A retailer
 * @author prussell
 *
 */
public interface Store {
	
	/**
	 * Get the name of the store
	 * @return Name of store
	 */
	public String getName();
	
	/**
	 * List the items sold at the store
	 */
	public void listItems();
	
	/**
	 * Get the price of an item
	 * @param item Item
	 * @return Price
	 */
	public float getPrice(Item item);
	
}
