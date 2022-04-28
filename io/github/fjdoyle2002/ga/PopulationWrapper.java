package io.github.fjdoyle2002.ga;

import io.github.fjdoyle2002.ga.threading.ScoringThreadPool;
import io.github.fjdoyle2002.util.RandomNumberGenerator;
import io.github.fjdoyle2002.util.io.SimpleLog;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;



public class PopulationWrapper {
	
	//private FitnessAlgorithm fitness_algorithm = null;
	private static SimpleLog log = null;
	private static SimpleLog census_file = null;
	private RandomNumberGenerator rndnumbergenerator = null;
	private Environment environment = null;	
	private ArrayList<Individual> population = null;
	private ArrayList<GenerationStatistics> statistics = null;
	private Set<Integer> persisted_individuals = null;
	private Map<Integer, Individual> individual_map = null;
	private int generation = 0;
	private ScoringThreadPool scoringThreadPool = null;
	
	/**
	 * Constructor
	 * @param environment
	 */
	public PopulationWrapper(Environment environment) {
		this.environment = environment;
		PopulationWrapper.log = Environment.getLog();
		this.statistics = new ArrayList<GenerationStatistics>();
		this.rndnumbergenerator = io.github.fjdoyle2002.util.RandomNumberGenerator.getInstance();
		this.persisted_individuals = new HashSet<Integer>();
		this.individual_map = new HashMap<Integer, Individual>();
	}
	
	/**
	 * Method runs the core evolutionary logic of the GA producing a successive generation 
	 */
	public void evolve(){
		this.generation++;
		Date now = new Date();
		log.comment("Beginning generation "+this.generation+" at "+now.toString()+".\n");
		double generational_mutation_severity = this.environment.getMaximum_mutation_severity()* rndnumbergenerator.getRandomDouble();		
		log.comment("Immigration ocurring...");
		this.immigrate();
		log.comment("Ranking current population...");		
		this.rank(this.population);
		this.sort();		
		log.comment("Generation mutation severity: "+generational_mutation_severity+"\n");
		log.comment("Introducing mutations to current population...");
		this.mutate(generational_mutation_severity);		
		log.comment("Breeding current population...");
		this.breed();		
		this.sort();
		this.take_census();
		log.comment("Population size prior to filtering: "+this.population.size());
		//filter prior to culling, to avoid removing unique members prior to redundant ones
		this.filterClones();
		log.comment("Population size post filtering: "+this.population.size());
		this.cull();		
		this.collect_stats();
		log.comment("Average scoring operation time in milliseconds: "+this.environment.getAverageScoreTime());
		this.logTopMembers();
		
	}
	private void logTopMembers() {
		log.comment("Generation "+this.generation+" top 10:");
		for(int i=1;i<=10;i++){
			log.comment("Rank ["+i+"]:"+this.population.get(this.population.size()-i));
		}
				
	}
	/**
	 * 
	 */
	public void cleanUp() {
		log.comment("Cleaning up threading if needed...");

		if(this.scoringThreadPool!=null){
			this.scoringThreadPool.shutDown();
			this.scoringThreadPool=null;
		}
		
	}

	/**
	 * 
	 */
	private void take_census() {
		Iterator<Individual> individual_iter = this.population.iterator();
		while(individual_iter.hasNext()){
			Individual current_individual = individual_iter.next();
			if(!this.persisted_individuals.contains(current_individual.getId())){
				recordIndividual(current_individual);
				//also, store in map for use in reporting
				individual_map.put(current_individual.getId(), current_individual);				
			}
		}
		
	}

	/**
	 * Method removes genetically redundant members of population
	 */
	private void filterClones() {
		Set<String>genotypes = new HashSet<String>();
		Iterator<Individual> pop_iter = this.population.iterator();
		while(pop_iter.hasNext()){
			Individual curr_individual = pop_iter.next();
			String curr_genotype = curr_individual.getStringRepresentation();
			if(genotypes.contains(curr_genotype)){
				//remove redundant individual
				pop_iter.remove();								
			}else{
				genotypes.add(curr_genotype);
			}			
		}
		
	}

	/**
	 * @return the current generation number
	 */
	public int getGeneration() {
		return generation;
	}
	
	/**
	 * Collects statistics about the current generation
	 */
	private void collect_stats() {
		GenerationStatistics gs = new GenerationStatistics(this.generation, this.population);
		log.comment("Stats: "+gs+"\n");
		this.statistics.add(gs);
		
	}
	/**
	 * Sorts current population by fitness
	 */
	private void sort() {
		Collections.sort(this.population);
		//Collections.reverse(this.population);
		
	}
	
	/**
	 * Removes least fit portion of population
	 */
	private void cull() {
		if(this.population.size()>this.environment.getMaximum_population_size()){
			ArrayList<Individual> subpopulation = new ArrayList<Individual>();
			int start_index = this.population.size() - this.environment.getMaximum_population_size();
			if(start_index > 0){			
				subpopulation.addAll(this.population.subList(start_index, (this.population.size())));
				this.population = subpopulation;
			}
		}
		
	}
	/**
	 * Method iterates over population members and calls scoring algorithm on each to allow sorting
	 * @param population
	 */
	private void rank(ArrayList<Individual> population) {
		
		if(this.environment.allowMultiThreading()){
			rankMultiThread(population);
		}else{
			
		
			Iterator<Individual> population_iterator = population.iterator();
			while(population_iterator.hasNext()){
				Individual current_individual = population_iterator.next();
				if(current_individual.needsScoring()){
					long start_time = System.currentTimeMillis();
					this.environment.getFitness_algorithm().score(current_individual, 1);
					long run_time = System.currentTimeMillis() - start_time;
					this.environment.registerScoreTime(run_time);
				}
			}
			
			
			
			
		}
		
		
		
	}
	/**
	 * 
	 * @param population
	 */
	private void rankMultiThread(ArrayList<Individual> population) {
		if(this.scoringThreadPool==null){
			this.scoringThreadPool=new ScoringThreadPool(this.environment);
		}
		List<Individual> individualsToBeScored = new ArrayList<Individual>();
		Iterator<Individual> population_iterator = population.iterator();
		while(population_iterator.hasNext()){
			Individual current_individual = population_iterator.next();
			if(current_individual.needsScoring()){
				individualsToBeScored.add(current_individual);
			}
		}
		
		this.scoringThreadPool.scoreIndividuals(individualsToBeScored);		
		
	}

	/**
	 * Introduces new group to population to prevent stagnation
	 */
	private void immigrate(){
		
		//List<Individual> foreign_pop = new ArrayList<Individual>();
		int migration_group_size = (int)Math.floor(environment.getImmigration_rate()*environment.getMaximum_population_size());
		for(int i=0;i<migration_group_size;i++){
			Individual new_member = environment.getIndividual_factory().getNewIndividual(this.generation);
			new_member.setOldest_ancestor(this.generation);
			new_member.setYoungest_ancestor(this.generation);
			new_member.setFullLineage("Rnd");
			new_member.setImmediateLineage("New Random Individual");
			this.population.add(new_member);						
		}
		
	}
	/**
	 * Recombines attributes of population members to produce offspring
	 */
	private void breed() {
		//determine how many breeding events to occur
		//determine who gets to breed high with medium, random? maybe percentages for h/h, m/h, m/m, l/h, l/m, l/l?
		int population_size = this.population.size();
		if(population_size>=2){
			ArrayList<Individual> children = new ArrayList<Individual>();
			double current_generational_fecundity=rndnumbergenerator.getRandomDouble()*this.environment.getGenerational_fecundity();
			int desired_offspring_count = (int)(Math.ceil((double)population_size * current_generational_fecundity)); 
			log.comment("Current generational fecundity: "+current_generational_fecundity+", desired_offspring_count: "+desired_offspring_count+"\n");
			while(children.size() < desired_offspring_count){
				int parent_1_index = (int)(Math.floor(rndnumbergenerator.getRandomDouble()*(double)population_size));			
				int parent_2_index = (int)(Math.floor(rndnumbergenerator.getRandomDouble()*(double)population_size));
				int oldest_ancestor = -1;
				int youngest_ancestor = -1;
				Individual parent_1 = population.get(parent_1_index);
				Individual parent_2 = population.get(parent_2_index);
				oldest_ancestor = parent_1.getOldest_ancestor() < parent_2.getOldest_ancestor() 
						? parent_1.getOldest_ancestor() 
							: parent_2.getOldest_ancestor();
				
				youngest_ancestor = parent_1.getGeneration() > parent_2.getGeneration()
						? parent_1.getGeneration()
							:parent_2.getGeneration();
			
				if(parent_1 != null && parent_2 != null && parent_1.getId() != parent_2.getId()){
					Individual[] offspring = parent_1.mate(parent_2);
					String full_lineage = "C[P1:"+parent_1.getId()+"|"+parent_1.getGeneration()+",L:"+parent_1.getFullLineage()+"]&[P2:G"+parent_2.getId()+"|"+parent_2.getGeneration()+",L:"+parent_2.getFullLineage()+"]";
					String immediate_lineage = "Child of["+parent_1.getId()+"&"+parent_2.getId()+"]";
					if(offspring[0] != null){
						offspring[0].setGeneration(this.generation);
						offspring[0].setOldest_ancestor(oldest_ancestor);
						offspring[0].setYoungest_ancestor(youngest_ancestor);
						offspring[0].setFullLineage(full_lineage);
						offspring[0].setImmediateLineage(immediate_lineage);
						children.add(offspring[0]);
						offspring[1].setGeneration(this.generation);
						offspring[1].setOldest_ancestor(oldest_ancestor);
						offspring[1].setYoungest_ancestor(youngest_ancestor);
						offspring[1].setFullLineage(full_lineage);
						offspring[1].setImmediateLineage(immediate_lineage);
						children.add(offspring[1]);												
					}
				}
			}
			log.comment("Ranking children...\n");
			this.rank(children);
			population.addAll(children);
					
		}
		
	}
	/**
	 * Introduces mutations to current population
	 * @param severity
	 */
	private void mutate(double severity) {
		Iterator<Individual> population_iterator = this.population.iterator();
		ArrayList<Individual> mutants = new ArrayList<Individual>();
		
		while(population_iterator.hasNext()){
			Individual current_individual = population_iterator.next();
			if(rndnumbergenerator.getRandomDouble() < severity){							
				Individual mutant = null;
				try {
					mutant = (Individual)current_individual.clone();					
					mutant.setGeneration(this.generation);
					//oldest ancestor should be set via clone
					mutant.mutate(severity);
					
					mutant.setFullLineage("M[id"+current_individual.getId()+"|"+current_individual.getGeneration()+",L:"+current_individual.getFullLineage()+"]");
					mutant.setImmediateLineage("Mutant of:["+current_individual.getId()+"]");
					mutants.add(mutant);

				} catch (CloneNotSupportedException e) {
					//this shouldn't happen, for now print stack trace and quit
					e.printStackTrace();
					System.exit(1);
				}				
			}
		}
		log.comment(mutants.size() + " mutants created. Ranking...\n");
		this.rank(mutants);
		
		this.population.addAll(mutants);
		
		
	}
	/**
	 * 
	 */
	public void setInitialPopulation(){
		this.population = new ArrayList<Individual>();
		for(int i=0;i<environment.getMaximum_population_size();i++){
			Individual new_member = environment.getIndividual_factory().getNewIndividual(this.generation);
			new_member.setOldest_ancestor(this.generation);
			new_member.setYoungest_ancestor(this.generation);
			new_member.setFullLineage("Rnd");
			new_member.setImmediateLineage("New Random Individual");
			this.population.add(new_member);						
		}
		log.comment("scoring initial population...");
		rank(this.population);
		
	}


	/**
	 * Returns the current population
	 * @return
	 */
	public ArrayList<Individual> getPopulation(){
		return this.population;
	}


	/**
	 * @return the statistics
	 */
	public ArrayList<GenerationStatistics> getStatistics() {
		return statistics;
	}


	/**
	 * 
	 * @param census
	 */
	public void setCensusFile(SimpleLog census){
		PopulationWrapper.census_file = census;
	}
	/**
	 * 
	 * @param individual
	 */
	private void recordIndividual(Individual individual){
		PopulationWrapper.census_file.comment(individual.getId()+"|FS:"+individual.getFitness_score()
				+"|Gen:"+individual.getGeneration()+"|OldestAncestor:"+individual.getOldest_ancestor()+"|YoungestAncestor:"+individual.getYoungest_ancestor()+"|AncOps:"+individual.getAncestral_operations()+"|Genes:"+individual.getStringRepresentation()
				+individual.getSupplementaryCensusData());
		this.persisted_individuals.add(individual.getId());
	}
	

}
