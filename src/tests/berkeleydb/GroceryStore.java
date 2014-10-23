package tests.berkeleydb;


import guttmanlab.core.util.StringParser;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;



import com.sleepycat.persist.model.Persistent;

/*
 * http://docs.oracle.com/cd/E17277_02/html/GettingStartedGuide/persistobject.html
 */
@Persistent
public class GroceryStore implements Store {
	
	private Map<Item, Float> prices;
	private String name;
	
	public GroceryStore() {}
	
	public GroceryStore(String storeName) {
		name = storeName;
		prices = new HashMap<Item, Float>();
	}
	
	/**
	 * Add an item to the store
	 * @param item Item
	 * @param price Price
	 */
	public void addItem(Item item, float price) {
		prices.put(item, Float.valueOf(price));
	}

	@Override
	public String getName() {
		return name;
	}

	@Override
	public void listItems() {
		System.out.println("\n" + name);
		for(Item item : prices.keySet()) {
			System.out.println(item.getName() + "\t" + item.getClass().getSimpleName() + "\t" + item.getColor() + "\t" + getPrice(item));
		}
	}

	/**
	 * Read from a table of stores and their inventory, and update items with the stores that sell them
	 * @param file Table in format: store_name   food_name   food_price
	 * @param itemsByName Items to update with the stores that sell them
	 * @throws IOException
	 */
	public static void readFromTable(String file, Map<String, ? extends Item> itemsByName) throws IOException {
		Map<String, GroceryStore> storesByName = new TreeMap<String, GroceryStore>();
		FileReader r1 = new FileReader(file);
		BufferedReader b1 = new BufferedReader(r1);
		StringParser s = new StringParser();
		while(b1.ready()) {
			s.parse(b1.readLine());
			if(s.getFieldCount() == 0) continue;
			String storeName = s.asString(0);
			@SuppressWarnings("unused")
			String foodName = s.asString(1);
			@SuppressWarnings("unused")
			float price = s.asFloat(2);
			if(!storesByName.containsKey(storeName)) {
				storesByName.put(storeName, new GroceryStore(storeName));
			}
			// Can't store instances of entity class (Food/Item)
			//storesByName.get(storeName).addItem(itemsByName.get(foodName), price);
			
		}
		r1.close();
		b1.close();
		FileReader r2 = new FileReader(file);
		BufferedReader b2 = new BufferedReader(r2);
		while(b2.ready()) {
			s.parse(b2.readLine());
			if(s.getFieldCount() == 0) continue;
			String storeName = s.asString(0);
			String foodName = s.asString(1);
			float price = s.asFloat(2);
			Store store = storesByName.get(storeName);
			itemsByName.get(foodName).addVendor(store, price);
		}
		r2.close();
		b2.close();
	}

	@Override
	public float getPrice(Item item) {
		return prices.get(item).floatValue();
	}
	
}
