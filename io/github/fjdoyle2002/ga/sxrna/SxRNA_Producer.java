package io.github.fjdoyle2002.ga.sxrna;
import io.github.fjdoyle2002.ga.GenerationStatistics;
import io.github.fjdoyle2002.ga.PopulationWrapper;
import io.github.fjdoyle2002.util.NonCodingRNA;
import io.github.fjdoyle2002.util.TargetRNA;
import io.github.fjdoyle2002.util.io.SequenceLoader;
import io.github.fjdoyle2002.util.io.SimpleLog;
import io.github.fjdoyle2002.util.io.TargetLoader;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Properties;






public class SxRNA_Producer {
	
	private SimpleLog log = null;
	
	private SxRNAEnvironment environment = null;
	private PopulationWrapper sxRNAPopulation = null;
	private String target_trigger = null;
	private String nuc_constraints = null;
	private String met_constraints = null;
	private String motif_struct = null;
	private String ncrna_filename = null;
	private String non_target_filename = null;
	private String secondary_positive_target_filename = null;
	private String negative_target_filename = null;
	private String properties_file = null;
	
	/**
	 * Constructor
	 */
	public SxRNA_Producer() {
		super();
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		SxRNA_Producer me = new SxRNA_Producer();
		me.init(args);
		me.runLogic();
	}
	
	
	private void runLogic() {
		int max_generation = this.environment.getMaximum_number_of_generations();
		for(int i=0;i<max_generation;i++){
		//for(int i=0;i<5;i++){
			this.sxRNAPopulation.evolve();
			//log.comment("\n-----------------Current Fittest----------------------------\n");
			//log.comment(this.sxRNAPopulation.getPopulation().get(this.sxRNAPopulation.getPopulation().size()-1));
		}
		this.sxRNAPopulation.cleanUp();
		//int highest_index = this.sxRNAPopulation.getPopulation().size()-1;
		log.comment("\n-----------------Stats----------------------------\n");
		List<GenerationStatistics> stats =this.sxRNAPopulation.getStatistics(); 
		Iterator<GenerationStatistics> stat_iter = stats.iterator();
		log.comment(GenerationStatistics.getHeader());
		while(stat_iter.hasNext()){
			GenerationStatistics currentStat = stat_iter.next();
			log.comment(currentStat.toString());
		}
	}
	/**
	 * 
	 * @param args
	 */
	private void init(String[] args) {
		this.parse(args);
		long time_stamp = System.currentTimeMillis();
		
		this.log = new SimpleLog("sxRNAProducer_"+time_stamp+".log");

		log.comment("initializing environment and population...");
		this.environment = new SxRNAEnvironment();
		SxRNAEnvironment.setLog(log);		
		if(null!=this.properties_file){
			loadPropertiesFromFile();
		}
		SxRNAEnvironment.setTarget_trigger(this.target_trigger);
		log.comment(this.environment.toString());
		
		processInputFiles();
	
		this.environment.setIndividual_factory(new SxRNAIndividualFactory());				
		
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
		
		this.sxRNAPopulation = new PopulationWrapper(this.environment);
		//this.sxRNAPopulation.setLog(log);
		SimpleLog census_file = new SimpleLog("sxRNACensus_"+time_stamp+".log");
		this.sxRNAPopulation.setCensusFile(census_file);
		//String nuc_constraints = "NNNNNNNNNNNNNNNNNNNNAAMGGYYCUUUUHAGRRCCMMNNNNNNNNNNNNNNNNNNNNN";
		//String met_constraints = "yyyyyyyyyyyyyyyyyyyynnnnnnnnnnnnnnnnnnnnnyyyyyyyyyyyyyyyyyyyyy";
		//String motif_struct =    ".......................((((((....)))))).......................";

		try {
			SxRNAIndividual.initializeStaticClassVariables(nuc_constraints, met_constraints, motif_struct);
		} catch (io.github.fjdoyle2002.util.sequence.BadFASTA_NucleotideSequenceException e) {
			log.comment("Bad FASTA exception initializing SxRNAIndividual static variables..."+e.getLocalizedMessage());
			System.exit(-1);
		}	
		//>hsa-miR-21-5p MIMAT0000076
		//UAGCUUAUCAGACUGAUGUUGA
		
		//mir122SxRNAIndividual.setTarget_trigger("UGGAGUGUGACAAUGGUGUUUG");
		//SxRNAIndividual.setTarget_trigger("UAGCUUAUCAGACUGAUGUUGA");
				
		

		
		//log.comment(non_targets);
		log.comment("setting initial population members...");
		this.sxRNAPopulation.setInitialPopulation();
		log.comment("SxRNAProducer initialized:"+this.toString());
				
	}
	/**
	 * Loads environment attributes for execution from a properties file
	 */
	private void loadPropertiesFromFile() {
		FileReader reader = null;
		try {
			reader = new FileReader(this.properties_file);
		} catch (FileNotFoundException e) {
			log.comment("Specified properties file does not exist: "+this.properties_file);
			System.exit(1);
		}  
	     Properties p=new Properties();  
		 try {
			p.load(reader);
		} catch (IOException e) {
			log.comment("IO error reading properties file: "+this.properties_file);
			System.exit(1);
			//e.printStackTrace();
		}
		this.environment.setAttributesFromProperties(p);
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
			}else if(args[i].startsWith("pf=")){
				this.properties_file = args[i].substring(3);
			}
			else{
				System.err.println("Unrecognized Commandline argument:"+args[i]);
				System.exit(1);
			}
		}
		
	}


	@Override
	public String toString() {
		return "SxRNA_Producer [log=" + log + ", motherNature=" + environment
				+ ", sxRNAPopulation=" + sxRNAPopulation + ", target_trigger="
				+ target_trigger + ", nuc_constraints=" + nuc_constraints
				+ ", met_constraints=" + met_constraints + ", motif_struct="
				+ motif_struct + ", ncrna_filename=" + ncrna_filename
				+ ", non_target_filename=" + non_target_filename
				+ ", secondary_positive_target_filename="
				+ secondary_positive_target_filename
				+ ", negative_target_filename=" + negative_target_filename
				+ "]";
	}



}
