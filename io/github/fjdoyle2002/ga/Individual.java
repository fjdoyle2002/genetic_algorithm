package io.github.fjdoyle2002.ga;

import io.github.fjdoyle2002.util.RandomNumberGenerator;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public abstract class Individual implements Cloneable, Comparable<Individual>{
	private static int idCount = 1;
	private int id =-1;
	private ArrayList<Gene> genes = null;
	private double fitness_score = 0;
	private RandomNumberGenerator rng = null;
	private boolean needsScoring = true;
	private static int[] default_crossover_points = null; //defines where reproductive recombination can take place by default
	private int generation = -1;
	private int oldest_ancestor = -1;
	private int youngest_ancestor = -1;
	private int ancestral_operations = 1; //init at 1, and increment at downstream operations
	private String stringRepresentation = null;
	private String full_lineage = null;
	private String immediate_lineage = null;
	private boolean notable_advance = false;
	private String note = "";
	
	/**
	 * Constructor
	 */
	protected Individual() {
		this.rng = io.github.fjdoyle2002.util.RandomNumberGenerator.getInstance();
		this.id = generateID();
	}
	/**
	 * 
	 * @return
	 */
	private static synchronized int generateID(){
		return Individual.idCount++;		
	}


	/**
	 * @return the fitness_score
	 */
	public double getFitness_score() {
		return fitness_score;
	}

	/**
	 * @param fitness_score the fitness_score to set
	 */
	public void setFitness_score(double fitness_score) {
		this.fitness_score = fitness_score;
		this.needsScoring = false;	
	}

	/**
	 * @param severity
	 */
	public void mutate(double severity){
		this.id = generateID();
		this.fitness_score = 0;
		this.setNeedsScoring(true);
		this.stringRepresentation = null;
		Iterator<Gene> gene_iterator = genes.iterator();
		while(gene_iterator.hasNext()){
			Gene current_gene = gene_iterator.next();
			if(this.rng.getRandomDouble() < severity){
				current_gene.mutate();				
			}
		}
		this.ancestral_operations++;
		this.postMutationHousekeeping();
		this.setStringRepresentation();
	}
	

	/**
	 * 
	 */
	protected abstract void postMutationHousekeeping();	
	/**
	 * 
	 * @param changedSinceScoring
	 */
	public void setNeedsScoring(boolean needsScoring) {
		this.needsScoring = needsScoring;
	}

	/**
	 * 
	 * @return
	 */
	public boolean needsScoring(){
		return this.needsScoring;		
	}

	/**
	 * Creates a child individual from combination of this and another parent Individual's genetic makeup
	 * @param parent_2
	 * @return
	 */
	public Individual[] mate(Individual parent_2){
		Individual[] offspring = new Individual[2];
		/*Must have legitimate mate*/
		if(parent_2 != null){
			int[] curr_crossover_points = this.getCrossover_points();
			/*Must have defined crossover points*/
			if(curr_crossover_points != null){
				int p2_ancestral_operations = parent_2.getAncestral_operations();
				int p1_ancestral_operations = this.getAncestral_operations();
				int offspring_ancestral_operations = p2_ancestral_operations + p1_ancestral_operations + 1;
				Individual child_1 = this.getChildInstance();
				child_1.setAncestral_operations(offspring_ancestral_operations);
				Individual child_2 = this.getChildInstance();
				child_2.setAncestral_operations(offspring_ancestral_operations);
				
				ArrayList<Gene> child_1_genes = new ArrayList<Gene>();
				ArrayList<Gene> child_2_genes = new ArrayList<Gene>();
				int last_locus = 0;
				int crossover_count = curr_crossover_points.length;
				for(int locus_index=0;locus_index <= crossover_count; locus_index++){
					int current_locus = 0;
					if(locus_index < crossover_count){
						current_locus = curr_crossover_points[locus_index];
					}else{
						current_locus = this.getGenes().size();						
					}
					double coin_toss = this.rng.getRandomDouble();
					//heads my genes to child 1
					List<Gene> genes_to_copy_1 = new ArrayList<Gene>(); 
					List<Gene> genes_to_copy_2 = new ArrayList<Gene>();
					if(coin_toss>.50){
						genes_to_copy_1.addAll(this.getGenes().subList(last_locus, current_locus));
						genes_to_copy_2.addAll(parent_2.getGenes().subList(last_locus, current_locus));
					}
					//tails my mate's genes to child 1
					else{
						genes_to_copy_2.addAll(this.getGenes().subList(last_locus, current_locus));
						genes_to_copy_1.addAll(parent_2.getGenes().subList(last_locus, current_locus));
					}
					Iterator<Gene> copy_iter = genes_to_copy_1.iterator();
					while(copy_iter.hasNext()){
						Gene curr_gene = copy_iter.next();
						try {
							child_1_genes.add((Gene)curr_gene.clone());
						} catch (CloneNotSupportedException e) {							
							e.printStackTrace();
							System.exit(1); //fatal error
						}
					}
					copy_iter = genes_to_copy_2.iterator();
					while(copy_iter.hasNext()){
						Gene curr_gene = copy_iter.next();
						try {
							child_2_genes.add((Gene)curr_gene.clone());
						} catch (CloneNotSupportedException e) {							
							e.printStackTrace();
							System.exit(1); //fatal error
						}
					}
					last_locus = current_locus;					
				}	
				child_1.setGenes(child_1_genes);			
				child_1.setStringRepresentation();
				child_2.setGenes(child_2_genes);			
				child_2.setStringRepresentation();
				
				offspring[0]=child_1;
				offspring[1]=child_2;
			}
			//else return array of nulls 
		}
		//else return array of nulls 
		return offspring;
	}
	
	
	
	/**
	 * 
	 * @return
	 */
	public abstract Individual getChildInstance();

	/**
	 * 
	 * @return
	 */
	public int getId(){
		return this.id;
	}

	/**
	 * 
	 * @return
	 */
	private  int[] getCrossover_points() {
		if(Environment.use_fixed_crossover()){
			return default_crossover_points;
		}else{
			int[] crossover_points = null;
			int xpoint1 = (int)(Math.ceil((double)this.genes.size() * this.rng.getRandomDouble()));  
			if(xpoint1 < 1){
				xpoint1 = 1;
			}
			
			double coin_toss = this.rng.getRandomDouble();
			if(coin_toss < .50){
				crossover_points = new int[1];				
				crossover_points[0] = xpoint1;
			}else{
				int xpoint2 = (int)(Math.ceil((double)this.genes.size() * this.rng.getRandomDouble()));
				while((Math.abs(xpoint2-xpoint1))<2){
					xpoint2 = (int)(Math.ceil((double)this.genes.size() * this.rng.getRandomDouble()));
				}
				crossover_points = new int[2];
				if(xpoint2<xpoint1){
					crossover_points[0]=xpoint2;
					crossover_points[1]=xpoint1;
				}else{
					crossover_points[0]=xpoint1;
					crossover_points[1]=xpoint2;
				}				
			}
			return crossover_points;
		}
	}
	/**
	 * 
	 * @return
	 */
	public ArrayList<Gene> getGenes() {
		return genes;
	}

	/**
	 * 
	 * @param genes
	 */
	protected void setGenes(ArrayList<Gene> genes) {
		this.genes = genes;
	}


	/**
	 * 
	 * @param crossover_points
	 */
	protected static void setDefaultCrossover_points(int[] crossover_points) {
		Individual.default_crossover_points = crossover_points;
	}


	/* (non-Javadoc)
	 * @see java.lang.Object#clone()
	 */
	@Override
	public Object clone() throws CloneNotSupportedException {
		Individual i2 = (Individual)super.clone();
		i2.genes = new ArrayList<Gene>();
		Iterator<Gene> geneIter = this.genes.iterator();
		while(geneIter.hasNext()){
			Gene currGene = geneIter.next();
			i2.genes.add((Gene)currGene.clone());			
		}
		return i2;
	}


	/* (non-Javadoc)
	 * @see java.lang.Comparable#compareTo(java.lang.Object)
	 */
	@Override
	public int compareTo(Individual o) {
		if(o == null){
			return 1;
		}else{
			if(o.fitness_score < this.fitness_score){
				return 1;
			}else if(o.fitness_score == this.fitness_score){
				return 0;
			}else{
				return -1;
			}
			
		}

	}


	/* (non-Javadoc)
	 * @see java.lang.Object#hashCode()
	 */
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		long temp;
		temp = Double.doubleToLongBits(fitness_score);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		return result;
	}


	/* (non-Javadoc)
	 * @see java.lang.Object#equals(java.lang.Object)
	 */
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		Individual other = (Individual) obj;
		if (Double.doubleToLongBits(fitness_score) != Double
				.doubleToLongBits(other.fitness_score))
			return false;
		return true;
	}

	public String getStringRepresentation() {
		if(this.stringRepresentation == null){
			this.setStringRepresentation();
		}
		return stringRepresentation;
	}

	/**
	 * Returns a string representation of the genetic makeup of this individual.
	 * No two individuals with distinct genetic makeup should return the same String from this
	 * method. This method is used to filter out genetic "clones" from population 
	 */
	public void setStringRepresentation() {
		StringBuffer sb = new StringBuffer();		
			if(this.genes != null){
				Iterator<Gene> geneIter = this.genes.iterator();
				while(geneIter.hasNext()){
					sb.append(geneIter.next());
				}
			}		
		this.stringRepresentation = sb.toString();
	}
	
	/**
	 * @return
	 */
	public int getGeneration() {
		return generation;
	}


	public void setGeneration(int generation) {
		this.generation = generation;
	}


	/**
	 * @return the full_lineage
	 */
	public String getFullLineage() {
		return full_lineage;
	}


	/**
	 * @param lineage the lineage to set
	 */
	public void setFullLineage(String full_lineage) {
		this.full_lineage = full_lineage;
	}

	/**
	 * 
	 * @return
	 */
	public String getImmediateLineage() {
		return immediate_lineage;
	}
	/**
	 * 
	 * @param immediate_lineage
	 */
	public void setImmediateLineage(String immediate_lineage) {
		this.immediate_lineage = immediate_lineage;
	}
	/**
	 * @return the notable_advance
	 */
	public boolean isNotable_advance() {
		return notable_advance;
	}


	/**
	 * @param notable_advance the notable_advance to set
	 */
	public void setNotable_advance(boolean notable_advance) {
		this.notable_advance = notable_advance;
	}


	/**
	 * @return the note
	 */
	public String getNote() {
		return note;
	}


	/**
	 * @param note the note to set
	 */
	public void setNote(String note) {
		this.note = note;
	}
	
	public abstract void setGenesByStringForTesting(String string_representation);
	/* (non-Javadoc)
	 * @see com.sxrnatechnologies.ga.Individual#toString()
	 */

	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		StringBuffer sb = new StringBuffer();
		sb.append("Individual ["+this.getId()+"], genes=" + this.genes + ", this.fitness_score="
				+ this.fitness_score);
		if(this.isNotable_advance()){
			sb.append(",***"+this.getNote()+"***");
		}
		sb.append(", needsScoring="
				+ this.needsScoring + ", generation=" + this.generation);
		sb.append(", Lineage:"+this.immediate_lineage);
		sb.append("\n");
		return sb.toString();
	}
	
	/**
	 * 
	 * @param oldest_ancestor
	 */
	public void setOldest_ancestor(int oldest_ancestor) {
		this.oldest_ancestor = oldest_ancestor;
	}
	/**
	 * 
	 * @return
	 */
	public int getOldest_ancestor() {
		return oldest_ancestor;
	}
	/**
	 * 
	 * @return
	 */
	public int getYoungest_ancestor() {
		return youngest_ancestor;
	}
	/**
	 * 
	 * @param youngest_ancestor
	 */
	public void setYoungest_ancestor(int youngest_ancestor) {
		this.youngest_ancestor = youngest_ancestor;
	}
	/**
	 * 
	 * @return
	 */
	public int getAncestral_operations() {
		return ancestral_operations;
	}
	/**
	 * 
	 * @param ancestral_operations
	 */
	public void setAncestral_operations(int ancestral_operations) {
		this.ancestral_operations = ancestral_operations;
	}	
	/**
	 * 
	 * @return
	 */
	public String getSupplementaryCensusData() {
		//by default, return nothing, if overriding class wants to 
		//add data, it is optional
		return "";
	}

	
}
