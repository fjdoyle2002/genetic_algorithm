package io.github.fjdoyle2002.ga.sxrna;

import io.github.fjdoyle2002.ga.Environment;
import io.github.fjdoyle2002.ga.FitnessAlgorithm;
import io.github.fjdoyle2002.ga.Individual;
import io.github.fjdoyle2002.parsers.RNAcofold_sxRNA_Parser;
import io.github.fjdoyle2002.parsers.SxRNA_Motif_Parser;
import io.github.fjdoyle2002.util.*;

import java.io.BufferedReader;
import java.io.File;

import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;




public class SxRNAFitnessAlgorithm implements FitnessAlgorithm {


	private static File working_directory = null;
	
	private static int NEGATIVE_TARGET_ENHANCEMENT = 1;
	private static int SECONDARY_POSITIVE_TARGET_ENHANCEMENT = 2;
	private static int NON_TARGET_ENHANCEMENT = 3;

	private static double sum_of_averages_for_non_targets = -1;
	private static double sum_of_fc_differentials_for_pos_targets = -1;
	private static double sum_of_fc_differentials_for_neg_targets = -1;
	
	private static double cleave_avoidance_bonus_factor = 10;
	private static double seed_binding_bonus_factor = 3;
	private static double t1_adenosine_avoidance_bonus = 8;
	
		 
	/**
	*/
	public SxRNAFitnessAlgorithm(SxRNAEnvironment environment) {	

		SxRNAFitnessAlgorithm.working_directory = new File(System.getProperty("user.dir"));
		SxRNAFitnessAlgorithm.initialize();

		
	}
	/**
	 * Perform any initial set up that the fitness algorithm requires
	*/
	private static void initialize() {

		setSumAveragesForNonTargets();
		setSumOfFCDifferentialsForPositiveTargets();
		setSumOfFCDifferentialsForNegativeTargets();
		
	}
	/**
	 * 
	 */
	private static void setSumOfFCDifferentialsForNegativeTargets() {
		Iterator<TargetRNA> neg_targ_iter =  SxRNAEnvironment.getNegative_targets().values().iterator();
		double sum = 0;
		while(neg_targ_iter.hasNext()){
			TargetRNA current = neg_targ_iter.next();
			double current_value = (1d/current.getOffFoldChange())*current.getNegTargetDifferentialValue();
			sum = sum + current_value;
		}
		SxRNAFitnessAlgorithm.sum_of_fc_differentials_for_neg_targets = sum;
		Environment.getLog().comment("sum_of_fc_differentials_for_neg_targets: "+sum_of_fc_differentials_for_neg_targets );
		
	}
	/**
	 * 
	 */
	private static void setSumOfFCDifferentialsForPositiveTargets() {
			
		Iterator<TargetRNA> pos_targ_iter =  SxRNAEnvironment.getSecondary_positive_targets().values().iterator();
		double sum = 0;
		while(pos_targ_iter.hasNext()){
			TargetRNA current = pos_targ_iter.next();
			double current_value = current.getOnFoldChange()*current.getPosTargetDifferentialValue();
			sum = sum + current_value;
		}
		SxRNAFitnessAlgorithm.sum_of_fc_differentials_for_pos_targets = sum;
		Environment.getLog().comment("sum_of_fc_differentials_for_pos_targets: "+sum_of_fc_differentials_for_pos_targets );
	}
	/**
	 * 
	 */
	private static void setSumAveragesForNonTargets() {
		Iterator<TargetRNA> non_targ_iter =  SxRNAEnvironment.getNon_targets().values().iterator();
		double sum = 0;
		while(non_targ_iter.hasNext()){
			TargetRNA current = non_targ_iter.next();
			sum = sum + current.getAverageValue();
			//System.out.println("Adding "+current.getAverageValue()+" to sum.");
		}
		SxRNAFitnessAlgorithm.sum_of_averages_for_non_targets = sum;
		Environment.getLog().comment("sum_of_averages_for_non_targets: "+sum_of_averages_for_non_targets );
	}

	@Override
	public double score(Individual individual, long thread_id) {
		double total_score = 0;
		SxRNAIndividual subject = null;
		try{
			subject = (SxRNAIndividual)individual;
		}catch (ClassCastException cce){
			System.err.println("Improper type sent to SxRNAFitness Algorithm...\n"+cce.getLocalizedMessage());
			
		}
		if(subject != null){
			double ncrna_interaction_scaling_factor = 0d;
			

			
			//check trigger independent motif presence
			List<String> fold_results = runFold(subject.getFASTA_representation());
			String fold = fold_results.get(1);
			subject.setRnaFold(fold);
			double percent_present = SxRNA_Motif_Parser.percentMotifPresenceInFold(fold, SxRNAIndividual.getMotifs(), subject.getBegin_motif_pos_in_fasta(), subject.getEnd_motif_pos_in_fasta(), SxRNAIndividual.getBegin_absolute_motif_struct_pos(),SxRNAIndividual.getEnd_absolute_motif_struct_pos());
			//calculate a penalty for this
			
			//use .5 threshold for declaring boolean of present in solo fold
			subject.setMotif_present_in_solo_fold(percent_present>=.5);
			
			
			double trigger_independent_score = 20 - (percent_present * 20);			
			double trigger_interaction_score = score_interaction_with_trigger(subject);
			
			
			double off_target_pos_interaction_score = 0;
			double off_target_neg_interaction_score = 0;
			if(SxRNAEnvironment.getNon_targets().size()>0){
				ncrna_interaction_scaling_factor = ncrna_interaction_scaling_factor+2;
				double[] results = score_interactions_with_other_triggers(SxRNAFitnessAlgorithm.NON_TARGET_ENHANCEMENT, subject, SxRNAEnvironment.getNon_targets().values());
				off_target_pos_interaction_score = results[0];
				off_target_neg_interaction_score = results[1];
				
			}
			
			
			double secondary_positive_target_interaction_score = 0;
			double secondary_positive_target_neg_interaction_score = 0;
			if(SxRNAEnvironment.getSecondary_positive_targets().size()>0){
				ncrna_interaction_scaling_factor = ncrna_interaction_scaling_factor+2;
				double[] results = score_interactions_with_other_triggers(SxRNAFitnessAlgorithm.SECONDARY_POSITIVE_TARGET_ENHANCEMENT, subject, SxRNAEnvironment.getSecondary_positive_targets().values());
				secondary_positive_target_interaction_score = results[0];
				secondary_positive_target_neg_interaction_score = results[1];
				
			}
			
			double negative_target_pos_interaction_score = 0;
			double targeted_neg_interaction_score = 0;
			if(SxRNAEnvironment.getNegative_targets().size()>0){
				ncrna_interaction_scaling_factor = ncrna_interaction_scaling_factor+2;
				double[] results = score_interactions_with_other_triggers(SxRNAFitnessAlgorithm.NEGATIVE_TARGET_ENHANCEMENT, subject, SxRNAEnvironment.getNegative_targets().values());
				negative_target_pos_interaction_score = results[0];
				targeted_neg_interaction_score = results[1];
				
			}
			if(ncrna_interaction_scaling_factor<1){
				ncrna_interaction_scaling_factor = 1;
			}
			total_score = (trigger_interaction_score * ncrna_interaction_scaling_factor) + (trigger_independent_score*ncrna_interaction_scaling_factor) + off_target_pos_interaction_score +secondary_positive_target_interaction_score
						+negative_target_pos_interaction_score + off_target_neg_interaction_score + targeted_neg_interaction_score + secondary_positive_target_neg_interaction_score;
			/*System.out.println("trigger_interaction_score: "+trigger_interaction_score);
			System.out.println("ncrna_interaction_scaling_factor: "+ncrna_interaction_scaling_factor);
			System.out.println("off_target_pos_interaction_score: "+off_target_pos_interaction_score);
			System.out.println("secondary_positive_target_interaction_score: "+secondary_positive_target_interaction_score);
			System.out.println("negative_target_pos_interaction_score: "+negative_target_pos_interaction_score);
			System.out.println("off_target_neg_interaction_score: "+off_target_neg_interaction_score);
			System.out.println("targeted_neg_interaction_score: "+targeted_neg_interaction_score);
			System.out.println("secondary_positive_target_neg_interaction_score: "+secondary_positive_target_neg_interaction_score);//*/
		}
		subject.setFitness_score(total_score);
		return total_score;
	}

	/**
	 * 
	 * @param subject
	 * @return
	 */
	private double[] score_interactions_with_other_triggers(int enhancement_type, SxRNAIndividual subject, Collection<TargetRNA> nonTargets) {
		double[] results = {0.0,0.0};
		//motif promoting interactions
		double pos_score = 0;
		double final_pos_score= 0;
		//suppressive interactions (not motif 'hiding', but message cleavage)
		double neg_score = 0;
		double final_neg_score= 0;
		
		//int pos_hit_count = 0;
		//int neg_hit_count = 0;
		
		StringBuffer pos_hits = new StringBuffer();
		StringBuffer neg_hits = new StringBuffer();
		
		if(nonTargets != null){
			Iterator<TargetRNA> nonTarget_iter = nonTargets.iterator();
			while(nonTarget_iter.hasNext()){
				
				TargetRNA current_nonTarget = nonTarget_iter.next();
				NonCodingRNA current_ncRNA = SxRNAEnvironment.retrieveNonCodingRNA(current_nonTarget);
				String cofold_input = null;
				
				cofold_input = subject.getFASTA_representation()+"&"+current_ncRNA.getSequence();
				
				//System.out.println("Checking cofold with "+cofold_input);
				List<String> fold_results = runCofold(cofold_input);
				double curr_neg_score = 0;
				if(fold_results.size()>1){
					curr_neg_score = checkNegativeInteraction(fold_results.get(1), subject, neg_hits, current_nonTarget, enhancement_type);
					if (curr_neg_score > 0){
						neg_score = neg_score + curr_neg_score;
						//neg_hit_count++;
					}
					//if an mirna shows suppressive potential, we don't bother checking for potential pos
					else{
						double curr_pos_score = checkPositiveInteraction(fold_results.get(1), subject, pos_hits, current_nonTarget, enhancement_type);						
						pos_score = pos_score + curr_pos_score;
						/*if(curr_pos_score>0){
							pos_hit_count++;
						}*/
					}
				}
				//try one more time before moving on. Do not treat as fatal error for one non-target fold, just log it
				else{
					System.err.println("Error checking non target fold for:"+cofold_input+", result size:"+fold_results.size());
					System.err.println("Making second attempt...");
					fold_results = runCofold(cofold_input);
					if(fold_results.size()>1){
						curr_neg_score = checkNegativeInteraction(fold_results.get(1), subject, neg_hits, current_nonTarget, enhancement_type);
						if (curr_neg_score > 0){
							neg_score = neg_score + curr_neg_score;
							//neg_hit_count++;
						}else{
							double curr_pos_score = checkPositiveInteraction(fold_results.get(1), subject, pos_hits, current_nonTarget, enhancement_type);						
							pos_score = pos_score + curr_pos_score;
							/*if(curr_pos_score>0){
								pos_hit_count++;
							}*/
						}						
					}else{
						System.err.println("Second error checking non target fold for:"+cofold_input+", result size:"+fold_results.size());
					}
				}
				
			}
		}
		//subject.setOffTargetPresence(hits.toString());
		if(enhancement_type==SxRNAFitnessAlgorithm.NEGATIVE_TARGET_ENHANCEMENT){			
			//score= SxRNAFitnessAlgorithm.sum_of_fc_differentials_for_neg_targets-score;			
			final_pos_score = 100-((pos_score/SxRNAFitnessAlgorithm.sum_of_fc_differentials_for_neg_targets)*100);
			subject.setNegTargetPosHits(pos_hits.toString());
			subject.setNegTargetPosHitScore(final_pos_score);
			final_neg_score = (neg_score/SxRNAFitnessAlgorithm.sum_of_fc_differentials_for_neg_targets)*100;
			subject.setNegTargetSuppressorHits(neg_hits.toString());
			subject.setNegTargetSuppressorHitScore(final_neg_score);
			
		}else if(enhancement_type==SxRNAFitnessAlgorithm.SECONDARY_POSITIVE_TARGET_ENHANCEMENT){
			final_pos_score = (pos_score/SxRNAFitnessAlgorithm.sum_of_fc_differentials_for_pos_targets)*100;
			subject.setPosTargetPosHits(pos_hits.toString());
			subject.setPosTargetPosHitScore(final_pos_score);
			final_neg_score = 100-((neg_score/SxRNAFitnessAlgorithm.sum_of_fc_differentials_for_pos_targets)*100);
			subject.setPosTargetSuppressorHits(neg_hits.toString());
			subject.setPosTargetSuppressorHitScore(final_neg_score);
					
		}else if(enhancement_type==SxRNAFitnessAlgorithm.NON_TARGET_ENHANCEMENT){
			final_pos_score = 100-((pos_score/SxRNAFitnessAlgorithm.sum_of_averages_for_non_targets)*100);
			subject.setNonTargetPosHits(pos_hits.toString());
			subject.setNonTargetPosHitScore(final_pos_score);
			final_neg_score = 100-((neg_score/SxRNAFitnessAlgorithm.sum_of_averages_for_non_targets)*100);
			subject.setNonTargetSuppressorHits(neg_hits.toString());
			subject.setNonTargetSuppressorHitScore(final_neg_score);
			
		}
		//System.out.println("In score_on_interactions... [score,final_score] reported as ["+score+","+final_score+"] for:"+ subject.getId());
		results[0] = final_pos_score;
		results[1] = final_neg_score;
		return results;
	}
	/**
	 * 
	 * @param cofold_result
	 * @param subject
	 * @param neg_hits
	 * @param current_nonTarget
	 * @param enhancement_type
	 * @return
	 */
	private double checkNegativeInteraction(String cofold_result, SxRNAIndividual subject, StringBuffer neg_hits, TargetRNA current_nonTarget, int enhancement_type) {
		double result = 0;
		NonCodingRNA current_ncRNA = SxRNAEnvironment.retrieveNonCodingRNA(current_nonTarget);
		boolean isSuppressing = false;
		if(SxRNAEnvironment.isTranslationalSwitch()){
			isSuppressing = RNAcofold_sxRNA_Parser.showsSuppressivePotential(cofold_result, subject.getFASTA_representation(), current_ncRNA.getSequence());
		}
		//if not a translational switch, just look to ensure the motif structure not substantially present in cofold
		else{
			
			isSuppressing = SxRNA_Motif_Parser.interferesWithMotif(cofold_result, subject.getBegin_motif_pos_in_fasta(), subject.getEnd_motif_pos_in_fasta());
		
		}
		
		
		if(isSuppressing){
			neg_hits.append(current_ncRNA.getId()+"|");
			if(enhancement_type==SxRNAFitnessAlgorithm.SECONDARY_POSITIVE_TARGET_ENHANCEMENT){
				result = current_nonTarget.getOnFoldChange()*current_nonTarget.getPosTargetDifferentialValue();
				neg_hits.append("["+current_nonTarget.getOn_state_value()+"]|");
			}else if(enhancement_type==SxRNAFitnessAlgorithm.NEGATIVE_TARGET_ENHANCEMENT){
				result = current_nonTarget.getOffFoldChange()*current_nonTarget.getNegTargetDifferentialValue();
				neg_hits.append("["+current_nonTarget.getOff_state_value()+"]|");
			}else if(enhancement_type==SxRNAFitnessAlgorithm.NON_TARGET_ENHANCEMENT){
				result = current_nonTarget.getAverageValue();
				neg_hits.append("["+current_nonTarget.getAverageScore()+"]|");
			}			
			
		}		
		return result;
	}
	/**
	 * 
	 * @param cofold_result
	 * @param subject
	 * @param pos_hits
	 * @param current_nonTarget
	 * @param enhancement_type
	 * @return
	 */
	private double checkPositiveInteraction(String cofold_result, SxRNAIndividual subject, StringBuffer pos_hits, TargetRNA current_nonTarget, int enhancement_type){
		double result = 0;
		NonCodingRNA current_ncRNA = SxRNAEnvironment.retrieveNonCodingRNA(current_nonTarget);
		boolean motif_present = RNAcofold_sxRNA_Parser.checkFoldForMotif(cofold_result, subject);
		if(motif_present){						
			pos_hits.append(current_ncRNA.getId()+"|");
			
			if(enhancement_type==SxRNAFitnessAlgorithm.SECONDARY_POSITIVE_TARGET_ENHANCEMENT){
				result = current_nonTarget.getOnFoldChange()*current_nonTarget.getPosTargetDifferentialValue();
				pos_hits.append("["+current_nonTarget.getOn_state_value()+"]|");
			}else if(enhancement_type==SxRNAFitnessAlgorithm.NEGATIVE_TARGET_ENHANCEMENT){
				result = (1d/current_nonTarget.getOffFoldChange())*current_nonTarget.getNegTargetDifferentialValue();
				pos_hits.append("["+current_nonTarget.getOff_state_value()+"]|");
			}else if(enhancement_type==SxRNAFitnessAlgorithm.NON_TARGET_ENHANCEMENT){
				result = current_nonTarget.getAverageScore();
				pos_hits.append("["+current_nonTarget.getAverageScore()+"]|");
			}
		}
		return result;
	}

	/**
	 * 
	 * @param subject
	 * @return
	 */
	private double score_interaction_with_trigger(SxRNAIndividual subject) {
		double trigger_interaction_score = 0;
		ThreeWayJunctionDescriptor results = checkForSplint(subject, SxRNAEnvironment.getTarget_trigger_name(), SxRNAEnvironment.getTarget_trigger());
		/*double binding_factor = results.getNcRNA_hybridization_percentage();*/	
		double spanning_factor = 0;
		double seed_bonus = 0;
		double cleavage_avoidance_bonus = 0;
		double dwell_avoidance_bonus = 0;
		if(results.ncrna_spans_motif()){
			spanning_factor = 1.0;
		}else{
			spanning_factor = 0.5;
		}
		//persist this portion of analysis in individual's properties 
		subject.setTrigger_spans_motif(results.ncrna_spans_motif());
		
		//double motif_factor = results.motif_presence();
		//trigger_interaction_score = (spanning_factor * (binding_factor * 60))+(motif_factor*40);
		double motif_factor = 0;
		if(results.motif_presence() == 1){
			motif_factor = 1;			
		}else{
			motif_factor = results.motif_presence() * .5;
		}

		
		//following are heuristic scoring bonuses based on literature for translational switches
		if(SxRNAEnvironment.isTranslationalSwitch()){
			int cleavage_multiplier = 3 - results.getNumber_of_critical_cleavage_bases_bound();
			if(cleavage_multiplier < 0){
				cleavage_multiplier = 0;
			}
			cleavage_avoidance_bonus = cleavage_multiplier * cleave_avoidance_bonus_factor;
			seed_bonus = results.getNumber_of_vital_seed_bases_bound() * seed_binding_bonus_factor;
			if(!results.isAdenosine_at_t1()){
				dwell_avoidance_bonus = t1_adenosine_avoidance_bonus;
			
			}			
		}
		//use 100% presence to set boolean of present in individual
		subject.setMotif_present_in_cofold(motif_factor==1);
		
		trigger_interaction_score = ((motif_factor * spanning_factor) * (results.getScore()+seed_bonus+cleavage_avoidance_bonus+dwell_avoidance_bonus));
		/*System.out.println("motif_factor:"+motif_factor);
		System.out.println("results.getScore():"+results.getScore());
		System.out.println("seed_bonus:"+seed_bonus);
		System.out.println("cleavage_avoidance_bonus:"+cleavage_avoidance_bonus);
		System.out.println("dwell_avoidance_bonus:"+dwell_avoidance_bonus);//*/
		
		//System.out.println("cofold score results: "+results);
		//System.out.println("trigger_interaction_score: "+trigger_interaction_score);
		return trigger_interaction_score;
	}
	
	
	public ThreeWayJunctionDescriptor checkForSplint(SxRNAIndividual candidate_RNA, String current_target_name, String current_target) {

		NonCodingRNA ncRNA = new NonCodingRNA(">"+current_target_name, current_target);
		ThreeWayJunctionDescriptor results = null;
		
		List<String> rawResults = runCofold(candidate_RNA.getFASTA_representation()+"&"+ncRNA.getSequence());
		  //now parse results and send to analysis module
		  //String cofoldResultString = parseCofoldResult();
		if(!(rawResults==null)&&rawResults.size()>1){
			candidate_RNA.setRnaCofold(rawResults.get(0)+"|"+rawResults.get(1));
			
			try{
				results = RNAcofold_sxRNA_Parser.evaluateThreeWayJunction(rawResults.get(1), candidate_RNA, ncRNA.getSequence());
				//System.out.println(rawResults);
			}catch(Exception e){
				System.err.println("Exception on call to evaluateThreeWayJunction()");
				System.err.println("Current combo: "+candidate_RNA +"&"+ncRNA.getSequence());
				System.err.println("RNAcofold result: "+rawResults.get(1));
				System.err.println("Motif Positions: "+(candidate_RNA.getBegin_motif_pos_in_fasta()) +", "+ (candidate_RNA.getEnd_motif_pos_in_fasta()));
				e.printStackTrace();
			}			
		}else{
			System.err.println("RNAcofold returned results of null or size 0 on full sequence query!");
			System.err.println(rawResults.get(1));
		}
		return results;
		
	}

	
	/**
	 *  Runs the cofold executable as defined by the "cofoldCommand" variable. 
	 */
	private List<String> runCofold(String sequences) {
		List<String> results = new ArrayList<String>();
		//Now execute PatSearch
		Process pr = null;
		BufferedReader input = null;
		BufferedReader error = null;
		try {
			
            Runtime rt = Runtime.getRuntime();
            
            //System.out.println("...calling: "+cofoldCommand);
            pr = rt.exec(SxRNAEnvironment.getCofoldCommand(), null, working_directory);
            PrintWriter pw = new PrintWriter(pr.getOutputStream());
            pw.println(sequences);
            pw.close();

            input =  new BufferedReader(new InputStreamReader(pr.getInputStream()));
            String standardOutput = input.readLine();
            while(standardOutput != null) {
            	results.add(standardOutput);
            	standardOutput=input.readLine();
            }
            error = new BufferedReader(new InputStreamReader(pr.getErrorStream()));
            String errorOutput =error.readLine();
            while(errorOutput != null) {
            	System.err.println(errorOutput);
                errorOutput=error.readLine();
            }  
           // System.out.println("...waiting for process end...");
            pr.waitFor();//don't continue until exit code received
            
            
        } catch(Exception e) {
            System.out.println(e.toString());
            e.printStackTrace();
        }
		//try to clean up the process resources	
        finally{
        	try {
				if (!(null==input)){input.close();}
			} catch (IOException e) {				
				e.printStackTrace();
			}
        	try {
				if (!(null==error)){error.close();}
			} catch (IOException e) {				
				e.printStackTrace();
			}	
			pr.destroy();
        }
        //System.out.println("Returning cofold results:"+results+" \nof size:"+results.size());
        return results;
        		
	}

	/**
	 *  Runs the fold executable as defined by the "foldCommand" variable. 
	 */
	public List<String> runFold(String sequence) {
		List<String> results = new ArrayList<String>();
		Process pr = null;
		BufferedReader input = null;
		BufferedReader error = null;
		try {
			//System.out.println("...getting runtime...");
            Runtime rt = Runtime.getRuntime();
            //Process pr = rt.exec(invocationString);
            //System.out.println("...calling: "+cofoldCommand);
            pr = rt.exec(SxRNAEnvironment.getFoldCommand(), null, working_directory);
            PrintWriter pw = new PrintWriter(pr.getOutputStream());
            pw.println(sequence);
            pw.close();

            input =  new BufferedReader(new InputStreamReader(pr.getInputStream()));
            String standardOutput = input.readLine();
            while(standardOutput != null) {
            	results.add(standardOutput);
            	standardOutput=input.readLine();
            }
            error = new BufferedReader(new InputStreamReader(pr.getErrorStream()));
            String errorOutput =error.readLine();
            while(errorOutput != null) {
            	System.err.println(errorOutput);
                errorOutput=error.readLine();
            }  
           // System.out.println("...waiting for process end...");
            pr.waitFor();//don't continue until exit code received
            
            
        } catch(Exception e) {
            System.out.println(e.toString());
            e.printStackTrace();
        }
			
        finally{
        	try {
				if (!(null==input)){input.close();}
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
        	try {
				if (!(null==error)){error.close();}
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}	
			pr.destroy();
        }
        //System.out.println("Returning cofold results:"+results+" \nof size:"+results.size());
        return results;
        		
	}	
	
	
}

