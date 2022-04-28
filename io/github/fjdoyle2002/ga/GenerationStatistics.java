 package io.github.fjdoyle2002.ga;

import java.util.ArrayList;
import java.util.Iterator;

public class GenerationStatistics {
	private int generation_count = 0;
	private int population_size = 0;
	private double _25th_percentile = 0;
	private double _50th_percentile = 0;
	private double _75th_percentile = 0;
	private double mean_score = 0;
	private double top_score = 0;

	public GenerationStatistics(int generation_count, ArrayList<Individual> population) {
		this.generation_count = generation_count;
		this.population_size = population.size();
		this.setStats(population);
	}
	/**
	 * Method takes a population of <Individual> objects and gathers basic statistics on score 
	 * @param population
	 */
	private void setStats(ArrayList<Individual> population){
		if(population!=null){
			Iterator<Individual> pop_iter = population.iterator();
			double score_sum = 0;
			double mean = 0;
			Individual current_individual = null;
			while(pop_iter.hasNext()){
				current_individual = pop_iter.next();
				score_sum = score_sum + current_individual.getFitness_score();				
			}
			this.top_score = current_individual.getFitness_score();
			
			mean = score_sum/(double)population.size();
			this.mean_score = mean;
			
			int _25th_percentile_pos = (int)(Math.round(.25*(population.size()+1)));
			this._25th_percentile = population.get(_25th_percentile_pos).getFitness_score();
			int _50th_percentile_pos = (int)(Math.round(.50*(population.size()+1)));
			this._50th_percentile = population.get(_50th_percentile_pos).getFitness_score();
			int _75th_percentile_pos = (int)(Math.round(.75*(population.size()+1)));
			this._75th_percentile = population.get(_75th_percentile_pos).getFitness_score();
			
		}
	}
	
		
	/**
	 * @return the generation_count
	 */
	public int getGeneration_count() {
		return generation_count;
	}
	/**
	 * @return the population_size
	 */
	public int getPopulation_size() {
		return population_size;
	}
	/**
	 * @return the _25th_percentile
	 */
	public double get_25th_percentile() {
		return _25th_percentile;
	}
	/**
	 * @return the _50th_percentile
	 */
	public double get_50th_percentile() {
		return _50th_percentile;
	}
	/**
	 * @return the _75th_percentile
	 */
	public double get_75th_percentile() {
		return _75th_percentile;
	}
	/**
	 * @return the mean_score
	 */
	public double getMean_score() {
		return mean_score;
	}
	/**
	 * @return the top_score
	 */
	public double getTop_score() {
		return top_score;
	}
	public static String getHeader() {
		return "Generation" + "\t "
				+ "population size" + "\t"
				+ "25th percentile" + "\t"
				+ "50th percentile" + "\t"
				+ "75th percentile" + "\t"
				+ "mean score" + "\t"
				+ "top score" + "\n";
	}
	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		return generation_count + "\t "
				+ population_size + "\t"
				+ _25th_percentile + "\t"
				+_50th_percentile + "\t"
				+ _75th_percentile+ "\t"
				+ mean_score + "\t"
				+ top_score + "\n";
	}


	
	
}

