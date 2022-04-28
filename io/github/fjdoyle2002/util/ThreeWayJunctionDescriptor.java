package io.github.fjdoyle2002.util;

import java.text.DecimalFormat;

/**
 * @author fd299212
 *
 */
public class ThreeWayJunctionDescriptor {
	private int lengthOfElongatedStem = 0;
	private double stem_elongation_percentage = 0;
	private double ncRNA_hybridization_percentage = 0;
	private double adjusted_ncRNAHybridization_percentage = 0;
	private double _5pBaitHelixBoundPercentage = 0;
	private double _3pBaitHelixBoundPercentage = 0;
	private double motif_presence = 0;
	private boolean ncrna_spans_motif = false;
	private int number_of_critical_cleavage_bases_bound = 0;
	private int number_of_vital_seed_bases_bound = 0;
	private boolean adenosine_at_t1 = false;
	//private double seed_region_binding_energy_estimate = 0.0;
	/**
	 * No args Constructor
	 */
	public ThreeWayJunctionDescriptor(){}
	/**
	 * Constructor
	 * @param ncRNAHybridizationPercentage
	 * @param adjusted_ncRNAHybridizationPercentage
	 * @param baitHelixBoundPercentage
	 * @param baitHelixBoundPercentage2
	 */
	public ThreeWayJunctionDescriptor(double fold_motif_presence, boolean ncrna_spans_motif, int lengthOfElongatedStem , double stem_elongation_percentage, double ncRNAHybridizationPercentage,
			double adjusted_ncRNAHybridizationPercentage,
			double _5pBaitHelixBoundPercentage, double _3pBaitHelixBoundPercentage, int vital_seed_bound, int critical_cleavage_bound, boolean adenosine_at_t1) {
		super();
		this.motif_presence = fold_motif_presence;
		this.ncrna_spans_motif = ncrna_spans_motif;
		this.stem_elongation_percentage = stem_elongation_percentage;
		this.ncRNA_hybridization_percentage = ncRNAHybridizationPercentage;
		this.adjusted_ncRNAHybridization_percentage = adjusted_ncRNAHybridizationPercentage;
		this._5pBaitHelixBoundPercentage = _5pBaitHelixBoundPercentage;
		this._3pBaitHelixBoundPercentage = _3pBaitHelixBoundPercentage;
		this.number_of_vital_seed_bases_bound = vital_seed_bound;
		this.number_of_critical_cleavage_bases_bound = critical_cleavage_bound;
		this.adenosine_at_t1 = adenosine_at_t1;
		//this.seed_region_binding_energy_estimate = mirna_seed_region_binding_energy;
	}
	/**
	 * @return the stem_elongation_percentage
	 */
	public double getStem_elongation_percentage() {
		return stem_elongation_percentage;
	}
	/**
	 * @param stem_elongation_percentage the stem_elongation_percentage to set
	 */
	public void setStem_elongation_percentage(double stem_elongation_percentage) {
		this.stem_elongation_percentage = stem_elongation_percentage;
	}	
	/**
	 * @return the ncRNAHybridizationPercentage
	 */
	protected double getNcRNAHybridizationPercentage() {
		return ncRNA_hybridization_percentage;
	}
	/**
	 * @return the adjusted_ncRNAHybridizationPercentage
	 */
	protected double getAdjusted_ncRNAHybridizationPercentage() {
		return adjusted_ncRNAHybridization_percentage;
	}
	/**
	 * @return the _5pBaitHelixBoundPercentage
	 */
	protected double get_5pBaitHelixBoundPercentage() {
		return _5pBaitHelixBoundPercentage;
	}
	/**
	 * @return the _3pBaitHelixBoundPercentage
	 */
	protected double get_3pBaitHelixBoundPercentage() {
		return _3pBaitHelixBoundPercentage;
	}
	
	public double getNcRNA_hybridization_percentage() {
		return ncRNA_hybridization_percentage;
	}
	public double getAdjusted_ncRNAHybridization_percentage() {
		return adjusted_ncRNAHybridization_percentage;
	}
	public double motif_presence() {
		return motif_presence;
	}
	public boolean ncrna_spans_motif() {
		return ncrna_spans_motif;
	}
	

	public double getScore(){
		double elongatedStemContribution = calculateElongatedStemContribution();
		double ncRNAHelixContribution = 100 - elongatedStemContribution;

		double elongateStemBasedScore = this.stem_elongation_percentage*elongatedStemContribution;
		double ncBasedScore = ((_5pBaitHelixBoundPercentage+_3pBaitHelixBoundPercentage)/2)*(adjusted_ncRNAHybridization_percentage)*ncRNAHelixContribution;
		double score = elongateStemBasedScore+ncBasedScore;
		return score;
	}
	/**
	 * 
	 * @return
	 */
	private double calculateElongatedStemContribution(){
		double result =0;
		if(this.lengthOfElongatedStem >0 && this.lengthOfElongatedStem < 10){			
			result = 10;
		}else if(this.lengthOfElongatedStem >10 && this.lengthOfElongatedStem < 20){
			result = 20;			
		}else if(this.lengthOfElongatedStem >20){
			result = 30;
		}
		return result;
	}
	/**
	 * 
	 * @return
	 */
	public int getNumber_of_critical_cleavage_bases_bound() {
		return number_of_critical_cleavage_bases_bound;
	}
	/**
	 * 
	 * @param number_of_vital_cleavage_bases_bound
	 */
	public void setNumber_of_critical_cleavage_bases_bound(
			int number_of_critical_cleavage_bases_bound) {
		this.number_of_critical_cleavage_bases_bound = number_of_critical_cleavage_bases_bound;
	}
	/**
	 * 
	 * @return
	 */
	public int getNumber_of_vital_seed_bases_bound() {
		return number_of_vital_seed_bases_bound;
	}
	/**
	 * 
	 * @param number_of_vital_seed_bases_bound
	 */
	public void setNumber_of_vital_seed_bases_bound(
			int number_of_vital_seed_bases_bound) {
		this.number_of_vital_seed_bases_bound = number_of_vital_seed_bases_bound;
	}
	/**
	 * 
	 * @return
	 */
	public boolean isAdenosine_at_t1() {
		return adenosine_at_t1;
	}
	/**
	 * 
	 * @param adenosine_at_t1
	 */
	public void setAdenosine_at_t1(boolean adenosine_at_t1) {
		this.adenosine_at_t1 = adenosine_at_t1;
	}
	/**
	 * 
	 * @param score
	 * @return
	 */
	private String formatDecimals(double score) {
		String scoreAsString = null;
		DecimalFormat scoreFormatter = new DecimalFormat("0.000");
		scoreAsString = scoreFormatter.format(score);
		return scoreAsString;
	}
	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		StringBuffer sb = new StringBuffer();
		sb.append("score:");
		sb.append(formatDecimals(this.getScore()));
		sb.append(" motif present: ");
		sb.append(this.motif_presence);
		sb.append(" 5p%:");
		sb.append(formatDecimals(this._5pBaitHelixBoundPercentage));
		sb.append(" 3p%:");
		sb.append(formatDecimals(this._3pBaitHelixBoundPercentage));
		sb.append(" nc%:");
		sb.append(formatDecimals(this.ncRNA_hybridization_percentage));
		sb.append(" adjNc%:");
		sb.append(formatDecimals(this.adjusted_ncRNAHybridization_percentage));
		sb.append(" stemElong%:");
		sb.append(formatDecimals(this.stem_elongation_percentage));
		sb.append(" number_of_vital_seed_bases_bound:");
		sb.append(this.number_of_vital_seed_bases_bound);
		sb.append(" number_of_vital_cleavage_bases_bound:");
		sb.append(this.number_of_critical_cleavage_bases_bound);
		sb.append(" adenosine at T1:");
		sb.append(this.adenosine_at_t1);
		
		return sb.toString();
	}

}
