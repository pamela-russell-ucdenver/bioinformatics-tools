package tests.berkeleydb;

import java.util.Iterator;

import broad.core.parser.CommandLineParser;

import com.sleepycat.persist.EntityCursor;

public class FoodStoreDelete {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-e", "Environment home directory", true);
		p.addStringArg("-f", "Food name to delete", false, null);
		p.addStringArg("-c", "Delete all foods with color", false, null);
		p.addBooleanArg("-a", "Delete all foods in database", false, false);
		p.parse(args);
		String envHome = p.getStringArg("-e");
		String foodName = p.getStringArg("-f");
		String color = p.getStringArg("-c");
		boolean all = p.getBooleanArg("-a");
		
		FoodDataAccessor a = new FoodDataAccessor(envHome, "test_store", false);
		
		if(foodName != null) {
			a.primaryIndex.delete(foodName);
		}
		
		if(color != null) {
			EntityCursor<Food> c = a.getAllFoodsWithColor(color);
			Iterator<Food> iter = c.iterator();
			while(iter.hasNext()) {
				Food f = iter.next();
				a.primaryIndex.delete(f.getName());
			}
			c.close();
		}
		
		if(all) {
			EntityCursor<Food> c = a.getAllFoods();
			Iterator<Food> iter = c.iterator();
			while(iter.hasNext()) {
				Food f = iter.next();
				a.primaryIndex.delete(f.getName());
			}
			c.close();
		}
		
		a.close();


	}

}
