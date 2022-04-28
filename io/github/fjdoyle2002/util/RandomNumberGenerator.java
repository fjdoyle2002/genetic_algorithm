package io.github.fjdoyle2002.util;

import ec.util.*;

public class RandomNumberGenerator {
	private static RandomNumberGenerator single_instance = null;
	
	private MersenneTwister mersenne_twister = null;
	
	/**
	 * Constructor
	 */
	private RandomNumberGenerator() {
		this.mersenne_twister = new MersenneTwister();
	
	}
	public synchronized static RandomNumberGenerator getInstance(){
		if(single_instance == null){
			single_instance = new RandomNumberGenerator();
		}
		return single_instance;
	}
	/**
	 * Returns a random double precision floating point between 0.0 and 1.0
	 * @return
	 */
	public synchronized double getRandomDouble(){
		return this.mersenne_twister.nextDouble();
	} 

	/**
	 * Replaces the instance of encapsulated random number generator with a new one.
	 */
	public synchronized void reset(){
		this.mersenne_twister = new MersenneTwister();
	}
}
