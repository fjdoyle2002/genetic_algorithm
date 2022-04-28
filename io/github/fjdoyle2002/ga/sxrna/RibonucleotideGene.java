package io.github.fjdoyle2002.ga.sxrna;

import io.github.fjdoyle2002.ga.Gene;
import io.github.fjdoyle2002.util.RandomNumberGenerator;


public class RibonucleotideGene extends Gene implements Cloneable{
	private boolean canBeBlank = false;
	private boolean mutable = false;
	private char[] possible_alleles = null;
	private char current_allele = ' ';
	private static final double[][] cutoffs = {{.00},{.50,1.0},{.33,.66,1.0},{.25,.50,.75,1.0},{.20,.40,.60,.80,1.0}};
	
	/**
	 * Constructor
	 * @param alleles
	 * @param canBeBlank
	 */
	public RibonucleotideGene(char[] alleles, boolean canBeBlank){
		super();
		this.canBeBlank = canBeBlank;
		this.possible_alleles = alleles;
		if(this.possible_alleles.length>1||this.canBeBlank){
			this.mutable = true;
		}else{
			this.current_allele = possible_alleles[0];
		}
		this.mutate();
		
	}
	/**
	 * 
	 * @param current_allele
	 */
	protected void setCurrentAlleleForTesting(char current_allele){
		this.current_allele = current_allele;
	}
	@Override
	public void mutate(){
		//only run logic if there are multiple choices for allele
		if(this.isMutable()){
			double random = RandomNumberGenerator.getInstance().getRandomDouble();
			//System.out.println("Mutate() rnd = "+random);
			int choices = this.possible_alleles.length-1;//adjust for indexing
			if(this.canBeBlank){
				choices=choices+1;//adjust for additional choice
			}
			double[] applicable_cutoffs = cutoffs[choices];
			int i=0;
			while(i<applicable_cutoffs.length){
				if(random <= applicable_cutoffs[i]){
					break;
				}					
				i++;
			}
			if(i>=this.possible_alleles.length){
				this.current_allele = ' ';
			}else{
				this.current_allele = this.possible_alleles[i];
			}
			//debug
			if((!this.canBeBlank) && this.current_allele == ' '){
				System.out.println("contradictory state:\n random = "+random+"\n current allele = |"+current_allele+"|\n i = "+i+"possible_alleles= "+this.possible_alleles);
			}
		}
	}
	/**
	 * 
	 * @return
	 */
	private final boolean isMutable(){
		return this.mutable;
	}
	
	/* (non-Javadoc)
	 * @see com.sxrnatechnologies.ga.Gene#clone()
	 */
	@Override
	protected Object clone() throws CloneNotSupportedException { 
		return super.clone();
	}

	/**
	 * 
	 */
	public String toString(){
		return String.valueOf(this.current_allele);
	}

}
