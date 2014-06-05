package tests.berkeleydb;

import java.io.File;

import nextgen.core.berkeleydb.DatabaseEnvironment;
import nextgen.core.berkeleydb.DatabaseStore;

import org.apache.log4j.Logger;


import com.sleepycat.je.DatabaseException;
import com.sleepycat.persist.EntityCursor;
import com.sleepycat.persist.PrimaryIndex;
import com.sleepycat.persist.SecondaryIndex;

/*
 * http://docs.oracle.com/cd/E17277_02/html/GettingStartedGuide/simpleda.html
 * http://docs.oracle.com/cd/E17277_02/html/GettingStartedGuide/persist_index.html
 */
/**
 * Data accessor which provides access to database
 * Opens a database environment
 * Must call close() when finished with this object
 * @author prussell
 *
 */
public class FoodDataAccessor {
	
	PrimaryIndex<String, Food> primaryIndex;
	SecondaryIndex<String, String, Food> secondaryIndexColor;
	private DatabaseEnvironment environment;
	private DatabaseStore entityStore;
	private static Logger logger = Logger.getLogger(FoodDataAccessor.class.getName());
	
	/**
	 * @param environmentHome Database environment home directory
	 * @param storeName Name of database store
	 * @param readOnly Whether environment should be read only
	 */
	public FoodDataAccessor(String environmentHome, String storeName, boolean readOnly) {
	
		/*
		 * http://docs.oracle.com/cd/E17277_02/html/GettingStartedGuide/persist_first.html
		 */
		environment = new DatabaseEnvironment();
		try {
			environment.setup(new File(environmentHome), readOnly, false);
		} catch(DatabaseException e) {
			logger.error("Problem with environment setup: " + e.toString());
			System.exit(-1);
		}
		
		entityStore = new DatabaseStore();
		try {
			entityStore.setup(environment.getEnvironment(), storeName, readOnly);
		} catch(DatabaseException e) {
			logger.error("Problem with store setup: " + e.toString());
			System.exit(-1);
		}
		
		primaryIndex = entityStore.getStore().getPrimaryIndex(String.class, Food.class);
		secondaryIndexColor = entityStore.getStore().getSecondaryIndex(primaryIndex, String.class, "color");

	}
	
	/**
	 * Close the store and the environment
	 */
	public void close() {
		entityStore.close();
		environment.close();
	}
	
	/**
	 * Get a cursor over all foods
	 * Must close when finished
	 * Can either use an iterator over the cursor or treat like a collection
	 * @return Cursor over all foods
	 */
	public EntityCursor<Food> getAllFoods() {
		return primaryIndex.entities();
	}
	
	/**
	 * Get a cursor over all foods with a particular color
	 * Must close when finished
	 * Can either use an iterator over the cursor or treat like a collection
	 * @param color Color
	 * @return Cursor over all foods with this color
	 */
	public EntityCursor<Food> getAllFoodsWithColor(String color) {
		return secondaryIndexColor.subIndex(color).entities();
	}
	
	/**
	 * Get a cursor over all foods with color in a particular range
	 * @param fromString String lower bound for color
	 * @param fromInclusive Lower bound is inclusive
	 * @param toString String upper bound for color
	 * @param toInclusive Upper bound is inclusive
	 * @return Cursor over all foods with colors in this interval
	 */
	public EntityCursor<Food> getAllFoodsWithColorBetween(String fromString, boolean fromInclusive, String toString, boolean toInclusive) {
		return secondaryIndexColor.entities(fromString, fromInclusive, toString, toInclusive);
	}
	
	
	/**
	 * Get a cursor over all foods with name in a particular range
	 * @param fromString String lower bound for name
	 * @param fromInclusive Lower bound is inclusive
	 * @param toString String upper bound for name
	 * @param toInclusive Upper bound is inclusive
	 * @return Cursor over all foods with names in this interval
	 */
	public EntityCursor<Food> getAllFoodsWithNameBetween(String fromString, boolean fromInclusive, String toString, boolean toInclusive) {
		return primaryIndex.entities(fromString, fromInclusive, toString, toInclusive);
	}
	
	
	
	/**
	 * Get a food by its primary key
	 * @param name Primary key value
	 * @return The food with this primary key
	 */
	public Food getFoodByName(String name) {
		return primaryIndex.get(name);
	}
	

	
}
