package tests.berkeleydb;

import guttmanlab.core.util.CommandLineParser;

import java.io.IOException;
import java.util.Collection;
import java.util.Map;

import org.apache.log4j.Logger;


/**
 * Class to add foods to the database
 * @author prussell
 *
 */
public class FoodStorePut {
	
	private FoodDataAccessor dataAccessor;
	@SuppressWarnings("unused")
	private static Logger logger = Logger.getLogger(FoodStorePut.class.getName());
	
	public FoodStorePut(String envHome, boolean readOnly) {
		dataAccessor = new FoodDataAccessor(envHome, "test_store", readOnly);
	}
	
	/*
	 * http://docs.oracle.com/cd/E17277_02/html/GettingStartedGuide/simpleput.html
	 */
	public void put(Collection<Food> foods) {
		for(Food food : foods) {
			dataAccessor.primaryIndex.put(food);
		}
	}
	
	public void put(Food food) {
		dataAccessor.primaryIndex.put(food);
	}
	
	/*
	 * http://docs.oracle.com/cd/E17277_02/html/GettingStartedGuide/simpleput.html
	 */
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-e", "Environment home directory", true);
		p.addBooleanArg("-r", "Database is read only", false, false);
		p.addStringArg("-f", "Table of foods (format: food_name   food_color)", true);
		p.addStringArg("-s", "Store table (format: store_name   food_name   food_price)", true);
		p.addStringArg("-i", "Table of ingredients (format: ingredient_name   ingredient_color   dish", true);
		p.parse(args);
		String envHome = p.getStringArg("-e");
		boolean readOnly = p.getBooleanArg("-r");
		String foodTable = p.getStringArg("-f");
		String storeTable = p.getStringArg("-s");
		
		FoodStorePut put = new FoodStorePut(envHome, readOnly);
		Map<String, ? extends Item> foods = Food.fromTable(foodTable);
		GroceryStore.readFromTable(storeTable, foods);
		for(Item item : foods.values()) {
			put.put((Food) item);
			Food food = put.dataAccessor.primaryIndex.get(item.getName());
			food.printSummary();
		}
				
		put.dataAccessor.close();
		
	}
	
}
