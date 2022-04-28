package io.github.fjdoyle2002.ga.orf;

import io.github.fjdoyle2002.ga.Gene;
import io.github.fjdoyle2002.util.RandomNumberGenerator;


public class CodonGene extends Gene{
	private boolean mutable = false;
	private String[] possible_alleles = null;
	private String current_allele = "";
	public CodonGene(String[] alleles){
		super();
		this.possible_alleles = alleles;
		if(this.possible_alleles.length>1){
			this.mutable = true;
			this.mutate();
		}else{
			this.current_allele = possible_alleles[0];
		}
		
		
	}
	
	@Override
	public void mutate(){
		//only run logic if there are multiple choices for allele
		if(this.isMutable()){
			double probability_segment = 1/(double)this.possible_alleles.length;	
			double random = RandomNumberGenerator.getInstance().getRandomDouble();
			double probability_sum = 0;
			int i;
			for(i=0;i<this.possible_alleles.length;i++){
				probability_sum = probability_sum+probability_segment;
				if(random < probability_sum){
					break;
				}
			}
			this.current_allele = this.possible_alleles[i];


		}
	}

	/**
	 * @return the mutable
	 */
	public boolean isMutable() {
		return mutable;
	}

	/**
	 * @param mutable the mutable to set
	 */
	public void setMutable(boolean mutable) {
		this.mutable = mutable;
	}

	/**
	 * @return the possible_alleles
	 */
	public String[] getPossible_alleles() {
		return possible_alleles;
	}

	/**
	 * @param possible_alleles the possible_alleles to set
	 */
	public void setPossible_alleles(String[] possible_alleles) {
		this.possible_alleles = possible_alleles;
	}

	/**
	 * @return the current_allele
	 */
	public String getCurrent_allele() {
		return current_allele;
	}

	/**
	 * @param current_allele the current_allele to set
	 */
	public void setCurrent_allele(String current_allele) {
		this.current_allele = current_allele;
	}
	/**
	 * 
	 */
	public String toString(){
		return this.current_allele;
	}
}
