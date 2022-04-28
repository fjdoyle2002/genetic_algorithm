package io.github.fjdoyle2002.parsers;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

public class SxRNA_Motif_Parser {

	
	/**
	 * Checks at what percentage any of the specified version of the motif structure are present in the fold.
     * Note, this only checks progressively smaller CENTERED portions of the motif*/
	 
	public static double percentMotifPresenceInFold(String foldOutput, List<String> MotifList, int specifiedSxRNAMotifStartIndex, int specifiedSxRNAMotifEndIndex, int motifDefinitionStartIndex, int motifDefinitionEndIndex){
		
		double presence = 0;
		//substring is beginning index inclusive, ending index exclusive...we want the char at the end index, so...
		int  corrected_sxrna_end_index = specifiedSxRNAMotifEndIndex+1;
		int corrected_motif_definition_end_index = motifDefinitionEndIndex+1;

		//because the evolved sxRNA sequence can vary in size due to "blank" nucleotide alleles, we have to 
		//differentiate between the motif structure position in the definition string and the motif related
		//ribonucleotide gene positions in the evolving sequence
		int adj_start_index = specifiedSxRNAMotifStartIndex;
		int adj_end_index = corrected_sxrna_end_index;
		int adj_motif_def_start_index = motifDefinitionStartIndex;
		int adj_motif_def_end_index = corrected_motif_definition_end_index;
		double full_motif_size = (double)corrected_sxrna_end_index - specifiedSxRNAMotifStartIndex;
		
		int partial_motif_size = adj_end_index - adj_start_index;
		
		outer_size_check_loop:
		while(partial_motif_size>0){
			
			double percent_of_total_size = (double)partial_motif_size/full_motif_size;
			String partial_motif_fold = foldOutput.substring(adj_start_index,adj_end_index);
			Iterator<String> motif_iter = MotifList.iterator();
			//we can't do a simple "contains" on the set as we are checking for position dependent partial match
			while(motif_iter.hasNext()){
				String current_motif = motif_iter.next();
				if(partial_motif_fold.equals(current_motif.subSequence(adj_motif_def_start_index, adj_motif_def_end_index))){
					presence = percent_of_total_size;
					break outer_size_check_loop;
				}
			}				
			adj_start_index++;
			adj_end_index--;
			adj_motif_def_start_index++;
			adj_motif_def_end_index--;
			partial_motif_size = adj_end_index - adj_start_index;	
			
		}
			
		
		
		return presence;
				
	}
	/**
	 * 
	 * @param foldOutput
	 * @param specifiedSxRNAMotifStartIndex
	 * @param specifiedSxRNAMotifEndIndex
	 * @return
	 */
	public static boolean interferesWithMotif(String foldOutput, int specifiedSxRNAMotifStartIndex, int specifiedSxRNAMotifEndIndex){
		boolean result = false;
		String parseLine = null;
		int energyPosition = foldOutput.indexOf("(-");
		//if we didn't see an energy reported in format'(-#'
		//check for it in '( -#' format
		if(energyPosition < 0){
			energyPosition = foldOutput.indexOf("( -");
		}		
		if(energyPosition > -1){
			parseLine = foldOutput.substring(0,energyPosition).trim();
		}
		else{
			parseLine = foldOutput.trim();
		}
		int[] positionMatches = DotParenthesesParser.parse(parseLine);
		int miRNAcutoffPosition = parseLine.indexOf('&');
		int interfering_bp_count = 0;
		int candidate_bp_count = 0;
		int nc_length = parseLine.length() - (miRNAcutoffPosition+1);
		for(int i=miRNAcutoffPosition+1;i< parseLine.length();i++){
			if(positionMatches[i] >=0 && positionMatches[i]<miRNAcutoffPosition){
				candidate_bp_count++;
				if(positionMatches[i] >= specifiedSxRNAMotifStartIndex && positionMatches[i]<=specifiedSxRNAMotifEndIndex){
					interfering_bp_count++;
				}				
			}
		}
		double percent_bound_to_candidate = (double)candidate_bp_count/(double)nc_length;
		//consider the interaction to be interfereing if more than 50% of the ncRNA is bound to candidateSxRNA
		//and any of the bp involve the motif region
		if(percent_bound_to_candidate > .5 && interfering_bp_count>0){
			result = true;
		}
		
		return result;
	}
	
	
}
