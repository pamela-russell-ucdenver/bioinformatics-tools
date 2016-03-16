package test.berkeleydb;

import guttmanlab.core.util.StringParser;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;


import com.sleepycat.persist.model.Entity;
import com.sleepycat.persist.model.PrimaryKey;
import com.sleepycat.persist.model.SecondaryKey;
import static com.sleepycat.persist.model.Relationship.*;

/*
 * http://docs.oracle.com/cd/E17277_02/html/GettingStartedGuide/persistobject.html
 */
@Entity(version=4)
public class Food implements Item {
	
	/*
	 * http://docs.oracle.com/cd/E17277_02/html/GettingStartedGuide/dplindexcreate.html
	 */
	@PrimaryKey
	private String name;
	@SecondaryKey(relate=MANY_TO_ONE)
	private String color;
	private Map<Store, Float> vendors;
	private Collection<String> dishes;
	
	public Food() {}
		
	public Food(String foodName, String foodColor, Collection<String> sampleDishes) {
		name = foodName;
		color = foodColor;
		vendors = new HashMap<Store, Float>();
		dishes = new TreeSet<String>();
		if(sampleDishes != null) {
			dishes.addAll(sampleDishes);
		}
	}
	
	public void addDish(String dish) {
		dishes.add(dish);
	}
	
	public void printSummary() {
		String summary = "\n";
		summary += "Ingredient name: " + name + "\n";
		summary += "Color: " + color + "\n";
		summary += "Vendors:\n";
		for(Store store : vendors.keySet()) {
			summary += store.getName() + " (" + vendors.get(store).toString() + ")\n";
		}
		summary += "Dishes:\n";
		for(String dish : dishes) {
			summary += dish + "\n";
		}
		System.out.println(summary);
	}
	
	public String getName() {
		return name;
	}
	
	public String getColor() {
		return color;
	}
	
	public float getPrice(Store store) {
		return store.getPrice(this);
	}
	
	public void setName(String n) {
		name = n;
	}
	
	public void setColor(String c) {
		color = c;
	}

	/**
	 * Line format: food_name   food_color   dish
	 * @param tableFile Table file
	 * @return Foods by name
	 * @throws IOException 
	 */
	public static Map<String, ? extends Item> fromTable(String tableFile) throws IOException {
		Map<String, Food> rtrn = new TreeMap<String, Food>();
		FileReader r = new FileReader(tableFile);
		BufferedReader b = new BufferedReader(r);
		StringParser s = new StringParser();
		while(b.ready()) {
			s.parse(b.readLine());
			if(s.getFieldCount() == 0) continue;
			String n = s.asString(0);
			String c = s.asString(1);
			String d = s.asString(2);
			if(!rtrn.containsKey(n)) {
				rtrn.put(n, new Food(n,c,null));
			}
			rtrn.get(n).addDish(d);
		}
		r.close();
		b.close();
		return rtrn;
	}

	@Override
	public void addVendor(Store store, float price) {
		vendors.put(store, Float.valueOf(price));
	}

	@Override
	public Map<Store, Float> getVendors() {
		return vendors;
	}
	
	
		
}
