package io.github.fjdoyle2002.ga;

public interface FitnessAlgorithm {

	public double score(Individual individual, long thread_id);
	
}
