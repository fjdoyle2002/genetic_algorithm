package io.github.fjdoyle2002.ga.sxrna;

import io.github.fjdoyle2002.ga.Individual;
import io.github.fjdoyle2002.ga.IndividualFactory;

public class SxRNAIndividualFactory implements IndividualFactory {

	@Override
	public Individual getNewIndividual(int generation) {
		SxRNAIndividual population_member = new SxRNAIndividual();
		population_member.setGeneration(generation);
		return population_member;
	}



}
