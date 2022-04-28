package io.github.fjdoyle2002.ga.sxrna;
import io.github.fjdoyle2002.ga.Individual;
import io.github.fjdoyle2002.util.NonCodingRNA;
import io.github.fjdoyle2002.util.TargetRNA;
import io.github.fjdoyle2002.util.io.SequenceLoader;
import io.github.fjdoyle2002.util.io.SimpleLog;
import io.github.fjdoyle2002.util.io.TargetLoader;

import java.io.File;
import java.util.Iterator;
import java.util.List;
import java.util.Map;


public class ScoreSxRNA {
	private SimpleLog log = null;
	
	private Individual individual_to_be_scored = null;
	private SxRNAEnvironment environment = null;
	private String target_trigger = null;
	private String nuc_constraints = null;
	private String met_constraints = null;
	private String motif_struct = null;
	private String ncrna_filename = null;
	private String non_target_filename = null;
	private String secondary_positive_target_filename = null;
	private String negative_target_filename = null;
	private String string_rep_individual_to_test = null;
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		ScoreSxRNA me = new ScoreSxRNA();
		me.init(args);
		me.runLogic();
	}
	
	
	private void runLogic() {
		long time_stamp = System.currentTimeMillis();
		log.comment("prior to scoring [at"+time_stamp+"...");
		log.comment(this.individual_to_be_scored);	
		this.environment.getFitness_algorithm().score(this.individual_to_be_scored, 1);
		time_stamp = System.currentTimeMillis();
		log.comment("after scoring [at"+time_stamp+"...");
		log.comment(this.individual_to_be_scored);
		
	}
	/**
	 * 
	 * @param args
	 */
	private void init(String[] args) {
		this.parse(args);
		long time_stamp = System.currentTimeMillis();
		
		this.log = new SimpleLog("sxRNAScorer_"+time_stamp+".log");
		
		log.comment("initializing environment and population...");
		this.environment = new SxRNAEnvironment();
		SxRNAEnvironment.setLog(log);
		SxRNAEnvironment.setTarget_trigger(this.target_trigger);
		
		processInputFiles();
		
		this.environment.setIndividual_factory(new SxRNAIndividualFactory());
		
		try {
			SxRNAIndividual.initializeStaticClassVariables(nuc_constraints, met_constraints, motif_struct);
		} catch (io.github.fjdoyle2002.util.sequence.BadFASTA_NucleotideSequenceException e) {
			log.comment("Bad FASTA exception initializing SxRNAIndividual static variables..."+e.getLocalizedMessage());
			System.exit(-1);
		}
						
		
		if(this.ncrna_filename !=null){
			this.environment.getEnvironmentalProperties().put("non_coding_rna_file", this.ncrna_filename);
		}
		if(this.non_target_filename !=null){
			this.environment.getEnvironmentalProperties().put("non_target_file", this.non_target_filename);
		}
		if(this.secondary_positive_target_filename !=null){
			this.environment.getEnvironmentalProperties().put("secondary_pos_target_file", this.secondary_positive_target_filename);
		}
		if( this.negative_target_filename!=null){
			this.environment.getEnvironmentalProperties().put("negative_target_file", this.negative_target_filename);			
		}
		
		this.environment.setFitness_algorithm(new SxRNAFitnessAlgorithm(environment));		
		
		
		//SxRNAEnvironment.setTarget_trigger(this.target_trigger);
		//this.environment.getProperties().put("non_target_file", this.non_target_filename);
		//this.environment.getProperties().put("secondary_pos_target_file", this.secondary_positive_target_filename);
		//this.environment.getProperties().put("negative_target_file", this.negative_target_filename);	
		
			
		log.comment("creating Individual to be scored");
		this.individual_to_be_scored = environment.getIndividual_factory().getNewIndividual(0);
		this.individual_to_be_scored.setGenesByStringForTesting(this.string_rep_individual_to_test);
		this.individual_to_be_scored.setStringRepresentation();
		log.comment("ScoreSxRNA initialized:"+this.toString());
		//*/
		
		
	}
	/**
	 * 
	 */
	private void processInputFiles() {
		File non_coding_rna_file = new File(this.ncrna_filename);
		Map<String, NonCodingRNA> nc_rna = SequenceLoader.loadNC_RNAs_as_map(non_coding_rna_file);
		log.comment("loaded "+nc_rna.size()+" non-coding RNAs...\n");
		SxRNAEnvironment.setNcRNAs(nc_rna);
		
		checkForTrigger(nc_rna);
		if(this.non_target_filename != null){
			File non_target_file = new File(this.non_target_filename);
			Map<String,TargetRNA> non_targets = TargetLoader.loadTarget_RNA_as_map(non_target_file);
			log.comment("loaded "+non_targets.size()+" non-targetRNAs...\n");
			SxRNAEnvironment.setNon_targets(non_targets);
		}
		
		if(this.secondary_positive_target_filename != null){
			File secondary_positive_target_file = new File(this.secondary_positive_target_filename);
			Map<String,TargetRNA> secondary_positive_targets = TargetLoader.loadTarget_RNA_as_map(secondary_positive_target_file);
			//if any of these is our specified trigger, remove as it is a special case in fitness calcs
			if(secondary_positive_targets.containsKey(SxRNAEnvironment.getTarget_trigger_name())){
				secondary_positive_targets.remove(SxRNAEnvironment.getTarget_trigger_name());
			}
			SxRNAEnvironment.setSecondary_positive_targets(secondary_positive_targets);
		}	
		if(this.negative_target_filename != null){
			File negative_target_file = new File(this.negative_target_filename);
			Map<String,TargetRNA> negative_targets = TargetLoader.loadTarget_RNA_as_map(negative_target_file);
			
			SxRNAEnvironment.setNegative_targets(negative_targets);
		}
		
	}
	/**
	 * 
	 * @param nc_rna
	 */
	private void checkForTrigger(Map<String, NonCodingRNA> nc_rna) {
		Iterator<NonCodingRNA> nc_iter = nc_rna.values().iterator();
		while(nc_iter.hasNext()){
			NonCodingRNA curr_nc_rna = nc_iter.next();
			if(SxRNAEnvironment.getTarget_trigger().equals(curr_nc_rna.getSequence())){
				SxRNAEnvironment.setTarget_trigger_name(curr_nc_rna.getId());
				break;
			}
		}
		
	}	
	/**
	 * Sets parameters from command line arguments
	 * @param args
	 */
	private void parse(String[] args) {
		for(int i=0;i<args.length;i++){
			if(args[i].startsWith("ms=")){
				this.motif_struct = args[i].substring(3);
			}else if(args[i].startsWith("mc=")){
				this.met_constraints = args[i].substring(3);				
			}else if(args[i].startsWith("nc=")){
				this.nuc_constraints = args[i].substring(3);				
			}else if(args[i].startsWith("tt=")){
				this.target_trigger = args[i].substring(3);				
			}else if(args[i].startsWith("ntf=")){
				this.non_target_filename= args[i].substring(4);
			}else if(args[i].startsWith("ncrf=")){
				this.ncrna_filename = args[i].substring(5);
			}else if(args[i].startsWith("sptf=")){
				this.secondary_positive_target_filename = args[i].substring(5);
			}else if(args[i].startsWith("netf=")){
				this.negative_target_filename = args[i].substring(5);

			}else if(args[i].startsWith("itt=")){
				this.string_rep_individual_to_test = args[i].substring(4);
			}else{
				System.err.println("Unrecognized Commandline argument:"+args[i]);
				System.exit(1);
			}
		}
		//System.out.println("ntf: "+this.non_target_filename);
		
	}





}
