package io.github.fjdoyle2002.ga;

import io.github.fjdoyle2002.util.io.SimpleLog;

import java.util.Properties;


public abstract class Environment {
	private static double maximum_mutation_severity = .5;
	private static double maximum_generational_fecundity = .5; 
	private static int maximum_population_size = 1000;
	private static double immigration_rate = .01;
	private static int maximum_number_of_generations = 100;
	private static FitnessAlgorithm fitness_algorithm = null;
	private static IndividualFactory individual_factory = null;
	private static boolean allowMultiThreading = true;
	private static long maximum_time_to_wait_for_thread = 25000;
	private static SimpleLog log = null;
	private static long score_times_sum = 0;
	private static long score_time_count = 0;	
	private static int maximum_threads = 32;
	private static boolean use_fixed_crossover = false;
	private static Properties environmental_properties = null;
	
	/**
	 * Constructor
	 */
	public Environment() {
		super();
		Environment.environmental_properties = new Properties();
	}
	/**
	 * Constructor
	 */
	public Environment(FitnessAlgorithm fa, IndividualFactory indfa) {
		super();
		Environment.fitness_algorithm = fa;
		Environment.individual_factory = indfa;
		Environment.environmental_properties = new Properties();
	}
	/**
	 * Accessor for the fitness algorithm 
	 * @return the fitness_algorithm set in this Environment object
	 * 		
	 */
	public FitnessAlgorithm getFitness_algorithm() {
		return fitness_algorithm;
	}
	/**
	 * Sets the fitness_algorithm for this Environment object
	 * @param fitness_algorithm the fitness_algorithm to set
	 */
	public void setFitness_algorithm(FitnessAlgorithm fitness_algorithm) {
		Environment.fitness_algorithm = fitness_algorithm;
	}
	/**
	 * Returns the immigration rate for the Environment object. This setting
	 * determines maximum number of new random candidates to introduce per generation 
	 * (i.e. per call to the evolve() method on PopulationWrapper class.
	 * @return the immigration_rate
	 */
	public double getImmigration_rate() {
		return immigration_rate;
	}
	/**
	 * Sets the immigration rate for the Environment object. This setting
	 * determines maximum number of new random candidates to introduce per generation 
	 * (i.e. per call to the evolve() method on PopulationWrapper class.
	 * @param immigration_rate the immigration_rate to set
	 */
	public void setImmigration_rate(double immigration_rate) {
		Environment.immigration_rate = immigration_rate;
	}
	/**
	 * 
	 * @return the maximum_number_of_generations
	 */
	public int getMaximum_number_of_generations() {
		return maximum_number_of_generations;
	}
	/**
	 * @param maximum_number_of_generations the maximum_number_of_generations to set
	 */
	public void setMaximum_number_of_generations(int maximum_number_of_generations) {
		Environment.maximum_number_of_generations = maximum_number_of_generations;
	}
	/**
	 * @return the maximum_population_size
	 */
	public int getMaximum_population_size() {
		return maximum_population_size;
	}
	/**
	 * @param maximum_population_size the maximum_population_size to set
	 */
	public void setMaximum_population_size(int maximum_population_size) {
		Environment.maximum_population_size = maximum_population_size;
	}
	/**
	 * @return the generational_fecundity
	 */
	public double getGenerational_fecundity() {
		return maximum_generational_fecundity;
	}
	/**
	 * @param generational_fecundity the generational_fecundity to set
	 */
	public void setGenerational_fecundity(double generational_fecundity) {
		Environment.maximum_generational_fecundity = generational_fecundity;
	}

	/**
	 * @return the maximum_mutation_severity
	 */
	public double getMaximum_mutation_severity() {
		return maximum_mutation_severity;
	}

	/**
	 * @param maximum_mutation_severity the maximum_mutation_severity to set
	 */
	public void setMaximum_mutation_severity(double maximum_mutation_severity) {
		Environment.maximum_mutation_severity = maximum_mutation_severity;
	}
	/**
	 * @return the individual_factory
	 */
	public IndividualFactory getIndividual_factory() {
		return individual_factory;
	}
	/**
	 * @param individual_factory the individual_factory to set
	 */
	public void setIndividual_factory(IndividualFactory individual_factory) {
		Environment.individual_factory = individual_factory;
	}
	public boolean allowMultiThreading() {
		return allowMultiThreading;
	}
	public void setAllowMultiThreading(boolean allowMultiThreading) {
		Environment.allowMultiThreading = allowMultiThreading;
	}
	public int getMaximum_threads() {
		return maximum_threads;
	}
	public void setMaximum_threads(int maximum_threads) {
		Environment.maximum_threads = maximum_threads;
	}
	public long getMaximum_time_to_wait_for_thread() {
		return maximum_time_to_wait_for_thread;
	}
	public void setMaximum_time_to_wait_for_thread(
			long maximum_time_to_wait_for_thread) {
		Environment.maximum_time_to_wait_for_thread = maximum_time_to_wait_for_thread;
	}
	public static SimpleLog getLog() {
		return Environment.log;
	}
	public static void setLog(SimpleLog log) {
		Environment.log = log;
	}
	public Properties getEnvironmentalProperties() {
		return environmental_properties;
	}
	public synchronized void registerScoreTime(long last_score_time) {
		Environment.score_times_sum = Environment.score_times_sum + last_score_time;
		Environment.score_time_count++;		
	}
	public double getAverageScoreTime(){
		double average_score_time = 0.0;
		average_score_time = (double)Environment.score_times_sum/(double)Environment.score_time_count;
		return average_score_time;
	}
	public static boolean use_fixed_crossover() {
		return use_fixed_crossover;
	}
	public static void setUse_fixed_crossover(boolean use_fixed_crossover) {
		Environment.use_fixed_crossover = use_fixed_crossover;
	}
	/**
	 * Uses values, if present, in properties file to override attribute default settings
	 */
	protected void setAttributesFromProperties(Properties properties) {
		Environment.environmental_properties.putAll(properties);
		if(Environment.environmental_properties.containsKey("maximum_mutation_severity")){
			try{
				double severity =  Double.parseDouble(Environment.environmental_properties.getProperty("maximum_mutation_severity"));
				if(severity >= 0.1 && severity <= .8){
					Environment.maximum_mutation_severity = severity;
				}else{
					Environment.log.comment("Invalid maximum_mutation_severity["+severity+"], using default..."+Environment.maximum_mutation_severity);
				}				
			}catch(NumberFormatException nfe){
				Environment.log.comment("Error parsing maximum_mutation_severity, using default..."+Environment.maximum_mutation_severity);
			}
		}
		if(Environment.environmental_properties.containsKey("maximum_generational_fecundity")){
			try{
				double fecundity = Double.parseDouble(Environment.environmental_properties.getProperty("maximum_generational_fecundity"));
				if(fecundity>=0.1 && fecundity <=.8){
					Environment.maximum_generational_fecundity = fecundity;
				}else{
					Environment.log.comment("Invalid maximum_generational_fecundity["+fecundity+"], using default..."+Environment.maximum_generational_fecundity);
				}
			}catch(NumberFormatException nfe){
				Environment.log.comment("Error parsing maximum_generational_fecundity, using default..."+Environment.maximum_generational_fecundity);
			}
		}
		if(Environment.environmental_properties.containsKey("immigration_rate")){
			try{
				double immigration = Double.parseDouble(Environment.environmental_properties.getProperty("immigration_rate"));
				if(immigration >= 0.0 && immigration <= .25){
					Environment.immigration_rate = immigration;
				}else{
					Environment.log.comment("Invalid immigration_rate["+immigration+"], using default..."+Environment.immigration_rate);
				}
			}catch(NumberFormatException nfe){
				Environment.log.comment("Error parsing immigration_rate, using default..."+Environment.immigration_rate);
			}
		}
		if(Environment.environmental_properties.containsKey("maximum_population_size")){
			try{
				int pop_size = Integer.parseInt(Environment.environmental_properties.getProperty("maximum_population_size"));
				if(pop_size >= 100 && pop_size <= 1000){
					Environment.maximum_population_size = pop_size; 
				}else{
					Environment.log.comment("Invalid maximum_population_size["+pop_size+"], using default..."+Environment.maximum_population_size);
				}
			}catch(NumberFormatException nfe){
				Environment.log.comment("Error parsing maximum_population_size, using default..."+Environment.maximum_population_size);
			}
		}
		if(Environment.environmental_properties.containsKey("maximum_number_of_generations")){
			try{
				int max_generations = Integer.parseInt(Environment.environmental_properties.getProperty("maximum_number_of_generations"));
				if(max_generations > 1 && max_generations < 100000){
					Environment.maximum_number_of_generations  = max_generations;
				}else{
					Environment.log.comment("Invalid maximum_number_of_generations["+max_generations+"], using default..."+Environment.maximum_number_of_generations);
				}
				
			}catch(NumberFormatException nfe){
				Environment.log.comment("Error parsing maximum_number_of_generations, using default..."+Environment.maximum_number_of_generations);
			}
		}
		if(Environment.environmental_properties.containsKey("maximum_threads")){
			try{
				Environment.maximum_threads = Integer.parseInt(Environment.environmental_properties.getProperty("maximum_threads"));
			}catch(NumberFormatException nfe){
				Environment.log.comment("Error parsing maximum_threads, using default..."+Environment.maximum_threads);
			}
		}
		if(Environment.environmental_properties.containsKey("maximum_time_to_wait_for_thread")){
			try{
				Environment.maximum_time_to_wait_for_thread = Long.parseLong(Environment.environmental_properties.getProperty("maximum_time_to_wait_for_thread"));
			}catch(NumberFormatException nfe){
				Environment.log.comment("Error parsing maximum_time_to_wait_for_thread, using default..."+Environment.maximum_time_to_wait_for_thread);
			}
		}
		if(Environment.environmental_properties.containsKey("allowMultiThreading")){
			try{
				Environment.allowMultiThreading = Boolean.parseBoolean(Environment.environmental_properties.getProperty("allowMultiThreading"));
			}catch(NumberFormatException nfe){
				Environment.log.comment("Error parsing allowMultiThreading, using default..."+Environment.allowMultiThreading);
			}
		}
		if(Environment.environmental_properties.containsKey("use_fixed_crossover")){
			try{
				Environment.use_fixed_crossover = Boolean.parseBoolean(Environment.environmental_properties.getProperty("use_fixed_crossover"));
			}catch(NumberFormatException nfe){
				Environment.log.comment("Error parsing use_fixed_crossover, using default..."+Environment.use_fixed_crossover);
			}
		}
		
	}
	@Override
	public String toString() {
		return "Environment [getImmigration_rate()=" + getImmigration_rate()
				+ ", getMaximum_number_of_generations()="
				+ getMaximum_number_of_generations()
				+ ", getMaximum_population_size()="
				+ getMaximum_population_size()
				+ ", getGenerational_fecundity()="
				+ getGenerational_fecundity()
				+ ", getMaximum_mutation_severity()="
				+ getMaximum_mutation_severity() + ", allowMultiThreading()="
				+ allowMultiThreading() + ", getMaximum_threads()="
				+ getMaximum_threads()
				+ ", getMaximum_time_to_wait_for_thread()="
				+ getMaximum_time_to_wait_for_thread() 
				+ ", use_fixed_crossover()=" + use_fixed_crossover()
				+"]";
	}

	

}
