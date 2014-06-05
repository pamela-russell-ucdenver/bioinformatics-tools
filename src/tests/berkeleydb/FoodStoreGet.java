package tests.berkeleydb;


import org.apache.log4j.Logger;

import com.sleepycat.persist.EntityCursor;

import broad.core.parser.CommandLineParser;

/**
 * Class to retrieve foods from the database by key values
 * @author prussell
 *
 */
public class FoodStoreGet {
	
	@SuppressWarnings("unused")
	private static Logger logger = Logger.getLogger(FoodStoreGet.class.getName());

	
	public static void main(String[] args) {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-e", "Environment home directory", true);
		p.addStringArg("-f", "Food name to get", false, null);
		p.addStringArg("-c", "Color to get", false, null);
		p.addBooleanArg("-a", "Print all foods in database", false, false);
		p.parse(args);
		String envHome = p.getStringArg("-e");
		String foodName = p.getStringArg("-f");
		String color = p.getStringArg("-c");
		boolean all = p.getBooleanArg("-a");
		
		FoodDataAccessor a = new FoodDataAccessor(envHome, "test_store", true);
		
		if(foodName != null) {
			Food food = a.getFoodByName(foodName);
			food.printSummary();
		}
		
		if(color != null) {
			EntityCursor<Food> c = a.getAllFoodsWithColor(color);
			for(Food f : c) {
				f.printSummary();
			}
			c.close();
		}
		
		if(all) {
			EntityCursor<Food> c = a.getAllFoods();
			for(Food f : c) {
				f.printSummary();
			}
			c.close();
		}
		
		a.close();
		
	}
	
}
