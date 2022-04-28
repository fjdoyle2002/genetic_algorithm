package io.github.fjdoyle2002.parsers;

import io.github.fjdoyle2002.ga.sxrna.SxRNAIndividual;
import io.github.fjdoyle2002.util.ThreeWayJunctionDescriptor;

import java.util.Iterator;



/**
 * @author fd299212
 * Class provides logic to analyze the trans three-way junction potential for a ncRNA vs an
 * RBP binding site via their RNAcofold results
 */
public class RNAcofold_sxRNA_Parser {
	/**
	 * Defines maximum allowed size of a three-way junction non bound junction region
	 * before penalty assessed 
	 */
	private static final int noPenaltyJunctionSize = 5;
		
	/**
	 * 
	 * @param cofoldOutput single string containing both mRNA and ncRNA fold, connected by '&'
	 * @param candidate_rna
	 * @param ncRNA_sequence 
	 * @return
	 */
	public static ThreeWayJunctionDescriptor evaluateThreeWayJunction(String cofoldOutput, SxRNAIndividual candidate_rna, String ncRNA_sequence) {

		    ThreeWayJunctionDescriptor twjd = null;
		//try{
			//first check if the required motif is present in the fold at position specified			
			
			//boolean foldContainsMotif = checkFoldForMotif(cofoldOutput, candidate_rna);
		    double motif_presence = SxRNA_Motif_Parser.percentMotifPresenceInFold(cofoldOutput, SxRNAIndividual.getMotifs(), candidate_rna.getBegin_motif_pos_in_fasta(), candidate_rna.getEnd_motif_pos_in_fasta(), SxRNAIndividual.getBegin_absolute_motif_struct_pos(),SxRNAIndividual.getEnd_absolute_motif_struct_pos());
			
			String parseLine = null;
			int energyPosition = cofoldOutput.indexOf("(-");
			//if we didn't see an energy reported in format'(-#'
			//check for it in '( -#' format
			if(energyPosition < 0){
				energyPosition = cofoldOutput.indexOf("( -");
			}
			
			if(energyPosition > -1){
				parseLine = cofoldOutput.substring(0,energyPosition).trim();
			}
			else{
				parseLine = cofoldOutput.trim();
			}
			
			int miRNAcutoffPosition = parseLine.indexOf('&');		 
			int _5pMotifBase = candidate_rna.getBegin_motif_pos_in_fasta();
			//int _5pStemBase = -1; //stem may continue past base of some motifs
			
			int _3pMotifBase = candidate_rna.getEnd_motif_pos_in_fasta();
			//int _3pStemBase = -1; //stem may continue past base of some motifs
			
			int lengthOfNcRNA = parseLine.length()-miRNAcutoffPosition-1;
			//variables to keep track of "bait" binding positions , initial values for ease of logic later
			//"first" refers to positions closest to motif 
			int firstBait5pSideComplexPos = 0;
			int lastBait5pSideComplexPos = 0;
			int firstBait3pSideComplexPos = 0;
			int lastBait3pSideComplexPos = 0;
			int ncRNA5pJunctionBasePos = 0;
			int ncRNA3pJunctionBasePos = 0;
			int ncRNAJunctionSize = 0;
			int _5pBaseOfElongatedStemPos = 0;
			int _3pBaseOfElongatedStemPos = 0;
			//double seed_region_binding_energy_estimate = 0;
			
			//variables for determining overall hybridization percentage of nucleotides in 
			//helices adjacent and near to motif
			double percentofBoundNucleotidesInElongatedStem = 0;
			int lengthOfElongateStem = 0;
			int numberOfBoundBaitNucleotidesIn5pHelix = 0;
			double percentOfBoundBaitNucleotidesIn5pHelix = 0;
			int numberOfBoundBaitNucleotidesIn3pHelix = 0;
			double percentOfBoundBaitNucleotidesIn3pHelix = 0;
			int[] positionMatches = DotParenthesesParser.parse(parseLine);
			//seed_region_binding_energy_estimate = calculateSeedEnergy(positionMatches, candidate_rna.getFASTA_representation(), ncRNA_sequence);
						
			boolean ncrna_spans_motif = ncRNAspansMotif(positionMatches, _5pMotifBase, _3pMotifBase, miRNAcutoffPosition);
			
				//first, get percentage of ncRNA that is bound to message
				int ncNucleotidesBoundToMessage = 0;
				for(int i=miRNAcutoffPosition+1;i<positionMatches.length;i++){
					//must be bound(>-1) to a nucleotide on the other strand
					int boundBasePos = positionMatches[i];
					if(boundBasePos>-1 && boundBasePos<miRNAcutoffPosition){
						ncNucleotidesBoundToMessage++;
						// we already know we are looking at a ncRNA base bound to a mRNA base
						if(boundBasePos > _3pMotifBase){
							numberOfBoundBaitNucleotidesIn3pHelix++;
							//set lastBait3pSideComplexPos to last mRNA nucleotide bound
							//if(lastBait3pSideComplexPos==0){
							if(boundBasePos > lastBait3pSideComplexPos){
								lastBait3pSideComplexPos = boundBasePos;
							}
							//set firstBait3pSideComplexPos to first mRNA nucleotide bound downstream of motif
							firstBait3pSideComplexPos = boundBasePos;
							ncRNA3pJunctionBasePos=i;
						}
						if(boundBasePos<_5pMotifBase){
							numberOfBoundBaitNucleotidesIn5pHelix++;
							if(boundBasePos > firstBait5pSideComplexPos){
								firstBait5pSideComplexPos=boundBasePos;
								ncRNA5pJunctionBasePos=i;
							}
							lastBait5pSideComplexPos = boundBasePos;
						}
					}
				}
				
				ncRNAJunctionSize = ncRNA5pJunctionBasePos-ncRNA3pJunctionBasePos-1;
				//System.out.println("ncRNA3pJunctionBasePos:"+ncRNA3pJunctionBasePos+" ncRNA5pJunctionBasePos:"+ncRNA5pJunctionBasePos+" ncRNAJunctionSize:"+ncRNAJunctionSize);
				double[] elongatedStemAnalysis = RNAcofold_sxRNA_Parser.analyzeElongatedBaitStem(positionMatches, firstBait5pSideComplexPos, firstBait3pSideComplexPos, _5pMotifBase, _3pMotifBase);
				_5pBaseOfElongatedStemPos = (int)elongatedStemAnalysis[0];
				_3pBaseOfElongatedStemPos = (int)elongatedStemAnalysis[1];
				lengthOfElongateStem = (int)elongatedStemAnalysis[2];
				percentofBoundNucleotidesInElongatedStem = elongatedStemAnalysis[3];
				
				int lengthOf5pBaitHelix = firstBait5pSideComplexPos-lastBait5pSideComplexPos+1;
				//int unbound_5p_junctionRegionToIncludeInPercentCalc = (_5pMotifBase-firstBait5pSideComplexPos-noPenaltyJunctionSize);
				int unbound_5p_junctionRegionToIncludeInPercentCalc = (_5pBaseOfElongatedStemPos-firstBait5pSideComplexPos-noPenaltyJunctionSize);
				if(unbound_5p_junctionRegionToIncludeInPercentCalc < 0){
					unbound_5p_junctionRegionToIncludeInPercentCalc=0;
				}
				percentOfBoundBaitNucleotidesIn5pHelix = (double)numberOfBoundBaitNucleotidesIn5pHelix/(double)(lengthOf5pBaitHelix+unbound_5p_junctionRegionToIncludeInPercentCalc);
				
				int lengthOf3pBaitHelix = lastBait3pSideComplexPos-firstBait3pSideComplexPos+1;
				//int unbound_3p_junctionRegionToIncludeInPercentCalc = (firstBait3pSideComplexPos-_3pMotifBase-noPenaltyJunctionSize);
				int unbound_3p_junctionRegionToIncludeInPercentCalc = (firstBait3pSideComplexPos-_3pBaseOfElongatedStemPos-noPenaltyJunctionSize);
				if(unbound_3p_junctionRegionToIncludeInPercentCalc < 0){
					unbound_3p_junctionRegionToIncludeInPercentCalc=0;
				}			
				percentOfBoundBaitNucleotidesIn3pHelix = (double)numberOfBoundBaitNucleotidesIn3pHelix/(double)(lengthOf3pBaitHelix+unbound_3p_junctionRegionToIncludeInPercentCalc);
				
				double percent_ncRNABoundToMessage=(double)ncNucleotidesBoundToMessage/(double)lengthOfNcRNA;
				int ncRNAJunctionSizeToRemoveFromPercentCalc = ncRNAJunctionSize;
				if(ncRNAJunctionSizeToRemoveFromPercentCalc>noPenaltyJunctionSize){
					ncRNAJunctionSizeToRemoveFromPercentCalc=noPenaltyJunctionSize;
				}
				double adjusted_percent_ncRNABoundToMessage=(double)ncNucleotidesBoundToMessage/(double)(lengthOfNcRNA-ncRNAJunctionSizeToRemoveFromPercentCalc);
				//now 	
				int vital_seed_bases_bound =  countVitalSeedBasesBound(positionMatches, candidate_rna.getFASTA_representation(), ncRNA_sequence);
				int critical_cleavage_bases_bound = countCriticalCleavagePositionsBound(positionMatches, candidate_rna.getFASTA_representation(), ncRNA_sequence);
				boolean adenosine_at_target_position_1 = checkTargetPositionOneForAdenosine(positionMatches, candidate_rna.getFASTA_representation(), ncRNA_sequence);
				
				twjd = new ThreeWayJunctionDescriptor(motif_presence, ncrna_spans_motif, lengthOfElongateStem, percentofBoundNucleotidesInElongatedStem, percent_ncRNABoundToMessage,adjusted_percent_ncRNABoundToMessage,percentOfBoundBaitNucleotidesIn5pHelix,percentOfBoundBaitNucleotidesIn3pHelix, vital_seed_bases_bound,critical_cleavage_bases_bound, adenosine_at_target_position_1) ;
			
			//}
		/*}catch(Exception e){
			System.err.println("Exception encountered in evaluateThreeWayJunction:"+e.getLocalizedMessage());
			System.err.println(cofoldOutput);
			System.err.println("specifiedMotifStartIndex: "+specifiedMotifStartIndex);
			System.err.println("specifiedMotifEndIndex: "+specifiedMotifEndIndex);
			e.printStackTrace();
			return null;
		}*/
		return twjd;
	}
	/**
	 * 
	 * @param cofoldOutput
	 * @param candidate_sequence
	 * @param ncRNA_sequence
	 * @return
	 */
	public static boolean showsSuppressivePotential(String cofoldOutput, String candidate_sequence, String ncRNA_sequence){
		boolean hasCleavagePotential = false;
		String parseLine = null;
		int energyPosition = cofoldOutput.indexOf("(-");
		//if we didn't see an energy reported in format'(-#'
		//check for it in '( -#' format
		if(energyPosition < 0){
			energyPosition = cofoldOutput.indexOf("( -");
		}
		
		if(energyPosition > -1){
			parseLine = cofoldOutput.substring(0,energyPosition).trim();
		}
		else{
			parseLine = cofoldOutput.trim();
		}				
		int[] positionMatches = DotParenthesesParser.parse(parseLine);
		
		if(countCriticalCleavagePositionsBound(positionMatches, candidate_sequence, ncRNA_sequence)> 1){
			if(calculateSeedEnergy(positionMatches, candidate_sequence, ncRNA_sequence) < -22 
							|| countVitalSeedBasesBound(positionMatches, candidate_sequence, ncRNA_sequence) >= 4){
				
				hasCleavagePotential = true;
			}
			
			
		}		
		return hasCleavagePotential;		
	}
	/**
	 * Estimates energy of base pairs for the miRNA "seed region"
	 * @param positionMatches
	 * @param fasta_representation
	 * @param ncRNA_sequence
	 * @return
	 */
	private static double calculateSeedEnergy(int[] positionMatches, String candidate_sequence, String ncRNA_sequence) {
		double energy = 0.0;
		int i = candidate_sequence.length();
		if(!(positionMatches[i]==-2)){
			System.err.println("Major parsing error in 'calculateSeedEnergy' for:"+candidate_sequence+":"+ncRNA_sequence);
			System.exit(-1);
		}
		String folded_sequence = candidate_sequence+"&"+ncRNA_sequence;
		int last_seed_position = i+7;
		if(positionMatches.length > last_seed_position){
			//adjust i to start at miRNA position 2
			for(i=i+2; i<=last_seed_position ;i++){
				if(positionMatches[i]>-1){
					energy = energy+getBasePairEnergy(folded_sequence.charAt(i),folded_sequence.charAt(positionMatches[i]));					
				}
			}		
		}
		return energy;
	}
	/**
	 * Determine number of bases on miRNA guide strand critical to cleavage are bound to target 
	 * @param positionMatches
	 * @param fasta_representation
	 * @param ncRNA_sequence
	 * @return
	 */
	private static int countCriticalCleavagePositionsBound(int[] positionMatches, String candidate_sequence, String ncRNA_sequence) {
		int positionsBound = 0;
		int i = candidate_sequence.length();
		if(!(positionMatches[i]==-2)){
			System.err.println("Major parsing error in 'calculateCriticalCleavagePositionsBound' for:"+candidate_sequence+":"+ncRNA_sequence);
			System.exit(-1);
		}
		String folded_sequence = candidate_sequence+"&"+ncRNA_sequence;
		int first_miRNA_base = i+1;
		int last_critical_position = i+11;
		if(positionMatches.length > last_critical_position){
			//adjust i to start at miRNA position 2
			for(i=i+9; i<=last_critical_position ;i++){
				//only count if bound and not to self
				if(positionMatches[i]>-1 && positionMatches[i]< first_miRNA_base ){
					positionsBound++;					
				}
			}		
		}
		return positionsBound;
	}	
	/**
	 * Determine if target position 1 on target strand is an "A" (affecting "dwell time"). This is done relative to first bound 
	 * base on guide strand after or including #2 (position 1 can't actually bind and what is predicted to be bound
	 * could have a bulge between it and next actual base pairing)  
	 * @param positionMatches
	 * @param fasta sequence of candidate
	 * @param ncRNA_sequence
	 * @return
	 */
	private static boolean checkTargetPositionOneForAdenosine(int[] positionMatches, String candidate_sequence, String ncRNA_sequence) {
		boolean hasAdenosineAtT1 = false;
		int i = candidate_sequence.length();
		if(!(positionMatches[i]==-2)){
			System.err.println("Major parsing error in 'checkTargetPositionOne' for:"+candidate_sequence+":"+ncRNA_sequence);
			System.exit(-1);
		}
		String folded_sequence = candidate_sequence+"&"+ncRNA_sequence;
		//position 'i' should be the ampersand connecting the two RNAs
		int mirna_strand_position_one = i+1;
		int first_mirna_paired_position = -1;
		//adjust i to start at miRNA position 2
		for(i=i+2; i<positionMatches.length ;i++){
			if(positionMatches[i]>-1){
				first_mirna_paired_position = i;
				break;					
			}
		}		
		if(first_mirna_paired_position >= 0){
			int paired_base = positionMatches[first_mirna_paired_position];
			//if the paired base is not a position less than the mirna start, the mirna is bound to itself
			if(paired_base < mirna_strand_position_one){
				int target_position_one = paired_base + (first_mirna_paired_position-mirna_strand_position_one);
				 char base_at_T1 = folded_sequence.charAt(target_position_one);
				 hasAdenosineAtT1 = (base_at_T1 == 'a' || base_at_T1 == 'A');
			}
		}
		
		return hasAdenosineAtT1;
	}


	/**
	 * Determine how many of guide strand seed positions 2 through 5 are bound to target as this has a substantial impact on
	 * how fast the miRNA finds it's target (should probably be checking if the bound target bases are reasonably 
	 * contiguous also...)  
	 * @param positionMatches
	 * @param fasta_representation
	 * @param ncRNA_sequence
	 * @return
	 */
	private static int countVitalSeedBasesBound(int[] positionMatches, String candidate_sequence, String ncRNA_sequence) {
		int positionsBound = 0;
		int i = candidate_sequence.length();
		if(!(positionMatches[i]==-2)){
			System.err.println("Major parsing error in 'calculateCriticalCleavagePositionsBound' for:"+candidate_sequence+":"+ncRNA_sequence);
			System.exit(-1);
		}
		String folded_sequence = candidate_sequence+"&"+ncRNA_sequence;
		int first_miRNA_base = i+1;
		int last_critical_position = i+5;
		if(positionMatches.length > last_critical_position){
			//adjust i to start at miRNA position 2
			for(i=i+2; i<=last_critical_position ;i++){
				//if bound and not bound to self
				if(positionMatches[i]>-1 && positionMatches[i]<first_miRNA_base){
					positionsBound++;					
				}
			}		
		}
		return positionsBound;
	}	
	/**
	 * Returns an energy value for a cis WC/WC bp of the two nucleotides provided
	 * @param base1
	 * @param base2
	 * @return
	 */
	private static double getBasePairEnergy(char base1, char base2) {
		double energy = 0.0;
		if(base1=='G'||base1=='g'){
			if(base2=='C'||base2=='c'){
				energy = -5.53;
			}else if(base2=='U'||base2=='u'){
				energy = -4.45;
			}
		}else if(base1=='C'||base1=='c'){
			if(base2=='G'||base2=='g'){
				energy = -5.53;
			}
		}else if(base1=='A'||base1=='a'){
			if(base2=='U'||base2=='u'){
				energy = -4.42;
			}			
		}else if(base1=='U'||base1=='u'){
			if(base2=='A'||base2=='a'){
				energy = -4.42;
			}else if(base2=='G'||base2=='g'){
				energy = -4.45;
			}
		}
		return energy;
	}
	/**
	 * Checks if the ncRNA portion of a RNACofold fold prediction spans the motif in the other portion
	 * @param positionMatches
	 * @param motifBase
	 * @param motifBase2
	 * @param miRNAcutoffPosition
	 * @return
	 */
	private static boolean ncRNAspansMotif(int[] positionMatches, int _5pMotifBase,
			int _3pMotifBase, int miRNAcutoffPosition) {
		
		//System.out.println("ncRNAspansMotif() entered with: "+ _5pMotifBase+" , "+_3pMotifBase+" miRNA_cutoof pos: "+miRNAcutoffPosition);
		boolean hits_5pSide = false;
		boolean hits_3pSide = false;
		//int i=0;
		int i = miRNAcutoffPosition+1;
		while(!(hits_5pSide && hits_3pSide) && i<positionMatches.length){
			if(positionMatches[i]>-1 && positionMatches[i]<_5pMotifBase){
				hits_5pSide = true;
			}
			if(positionMatches[i]>-1&&positionMatches[i]>_3pMotifBase){
				hits_3pSide = true;
			}
			i++;
		}
		return (hits_5pSide&&hits_3pSide);
	}
	/**
	 * 
	 * @param positionMatches
	 * @param _5pJunctionPos
	 * @param _3pJunctionPos
	 * @param _5pMotifBasePos
	 * @param _3pMotifBasePos
	 * @return 
	 */
	private static double[] analyzeElongatedBaitStem(int[] positionMatches, int _5pJunctionPos, int _3pJunctionPos, int _5pMotifBasePos, int _3pMotifBasePos){
		int _5pElongatedStemBase = _5pMotifBasePos; //init to base of motif
		int _3pElongatedStemBase = _3pMotifBasePos; //init  to base of motif
		double percentElongatedStemBound = 0;
		int totalElongatedStemSize = 0;
		/* first, we need to find the first bound nucleotide after junction regions prior to motif on each side*/
		for(int i=_5pJunctionPos+1;i<_5pMotifBasePos;i++){
			if(positionMatches[i] > 0){
				if(i < _5pElongatedStemBase){
					_5pElongatedStemBase = i;
				}
			}
		}
		for(int j=_3pMotifBasePos+1;j<_3pJunctionPos;j++){
			if(positionMatches[j] > 0){
				if(j > _3pElongatedStemBase){
					_3pElongatedStemBase = j;
				}
			}			
		}
		//if there is an elongated stem between motif and junction region, we need to address it
		if(_5pElongatedStemBase < _5pMotifBasePos && _3pElongatedStemBase > _3pMotifBasePos){
			/* Then, we need to determine how many nucleotides between those just identified and the base of the motif are bound*/
			int boundCount = 0;			
			for(int k=_5pElongatedStemBase;k<_5pMotifBasePos;k++){
				totalElongatedStemSize++;
				if(positionMatches[k] > 0){
					boundCount++;
				}
			}
			for(int l=_3pMotifBasePos+1;l<=_3pElongatedStemBase;l++){
				totalElongatedStemSize++;
				if(positionMatches[l] > 0){
					boundCount++;
				}
			}
			if(totalElongatedStemSize>0){
				percentElongatedStemBound = (double)boundCount/(double)totalElongatedStemSize;
			}
		}
		return new double[]{(double)_5pElongatedStemBase, (double)_3pElongatedStemBase, (double)totalElongatedStemSize, percentElongatedStemBound };
	}
	/**
	 * Checks presence of desired motif in fold
	 * @param foldOutput
	 * @param candidate
	 * @return
	 */
	public static boolean checkFoldForMotif(String foldOutput, SxRNAIndividual candidate) {
		boolean present = false;
		//substring is beginning index inclusive, ending index exclusive...we want the char at the end index, so...
		int  corrected_sxrna_end_index = candidate.getEnd_motif_pos_in_fasta()+1;
		int corrected_motif_definition_end_index = SxRNAIndividual.getEnd_absolute_motif_struct_pos()+1;

		//because the evolved sxRNA sequence can vary in size due to "blank" nucleotide alleles, we have to 
		//differentiate between the motif structure position in the definition string and the motif related
		//ribonucleotide gene positions in the evolving sequence
		int adj_start_index = candidate.getBegin_motif_pos_in_fasta();
		int adj_end_index = corrected_sxrna_end_index;
		int adj_motif_def_start_index = SxRNAIndividual.getBegin_absolute_motif_struct_pos();
		int adj_motif_def_end_index = corrected_motif_definition_end_index;
			
	
		String partial_motif_fold = foldOutput.substring(adj_start_index,adj_end_index);
		Iterator<String> motif_iter = SxRNAIndividual.getMotifs().iterator();
		//we can't do a simple "contains" on the set as we are checking for position dependent partial match
		while(motif_iter.hasNext()){
			String current_motif = motif_iter.next();
			if(partial_motif_fold.equals(current_motif.subSequence(adj_motif_def_start_index, adj_motif_def_end_index))){
				present = true;
				break;
			}
		}					
		
		return present;
	}	
	
	
}
