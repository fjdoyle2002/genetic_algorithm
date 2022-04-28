package io.github.fjdoyle2002.ga.sxrna;

import io.github.fjdoyle2002.ga.Environment;
import io.github.fjdoyle2002.util.NonCodingRNA;
import io.github.fjdoyle2002.util.TargetRNA;

import java.util.HashMap;
import java.util.Map;
import java.util.Properties;


/**
 * @author fd299212
 *
 */
public class SxRNAEnvironment extends Environment {
	private static String target_trigger = null;
	private static String target_trigger_name = "target_trigger";
	private static Map<String,NonCodingRNA> ncRNAs = null;
	private static Map<String,TargetRNA> non_targets = null;
	private static Map<String,TargetRNA> secondary_positive_targets = null;
	private static Map<String,TargetRNA> negative_targets = null;
	
	private static String cofoldCommand = "/usr/local/bin/RNAcofold -d2 --noLP --noPS";
	private static String foldCommand = "/usr/local/bin/RNAfold -p -d2 --noLP --noPS --noDP";
			
	
	private static boolean isTranslationalSwitch = true;
	
	public SxRNAEnvironment(){
		super();
		//initialize these...they will likely be replaced rather than added to, but 
		//easier for now than adding conditional logic where they are used FJD 7/2021
		SxRNAEnvironment.non_targets = new HashMap<String,TargetRNA>();
		SxRNAEnvironment.secondary_positive_targets = new HashMap<String,TargetRNA>();
		SxRNAEnvironment.negative_targets = new HashMap<String,TargetRNA>();
	}
	
	public static String getTarget_trigger() {
		return target_trigger;
	}
	public static void setTarget_trigger(String target_trigger) {
		SxRNAEnvironment.target_trigger = target_trigger;
	}
	
	public static String getTarget_trigger_name() {
		return target_trigger_name;
	}


	public static void setTarget_trigger_name(String target_trigger_name) {
		SxRNAEnvironment.target_trigger_name = target_trigger_name;
	}


	public static Map<String, NonCodingRNA> getNcRNAs() {
		return ncRNAs;
	}
	public static void setNcRNAs(Map<String, NonCodingRNA> ncRNAs) {
		SxRNAEnvironment.ncRNAs = ncRNAs;
	}

	public static Map<String, TargetRNA> getNon_targets() {
		return non_targets;
	}


	public static void setNon_targets(Map<String, TargetRNA> non_targets) {
		SxRNAEnvironment.non_targets = non_targets;
	}


	public static Map<String, TargetRNA> getSecondary_positive_targets() {
		return secondary_positive_targets;
	}


	public static void setSecondary_positive_targets(
			Map<String, TargetRNA> secondary_positive_targets) {
		SxRNAEnvironment.secondary_positive_targets = secondary_positive_targets;
	}


	public static Map<String, TargetRNA> getNegative_targets() {
		return negative_targets;
	}


	public static void setNegative_targets(Map<String, TargetRNA> negative_targets) {
		SxRNAEnvironment.negative_targets = negative_targets;
	}

	protected static NonCodingRNA retrieveNonCodingRNA(TargetRNA target){
		NonCodingRNA ncRNA = null;
		ncRNA = SxRNAEnvironment.ncRNAs.get(target.getId());
		return ncRNA;
	}

	public static String getCofoldCommand() {
		return SxRNAEnvironment.cofoldCommand;
	}
	public static String getFoldCommand() {
		return SxRNAEnvironment.foldCommand;
	}
	public static boolean isTranslationalSwitch() {		
		return SxRNAEnvironment.isTranslationalSwitch ;
	}

	public void setAttributesFromProperties(Properties properties){		
		super.setAttributesFromProperties(properties);
		if(properties.containsKey("isTranslationalSwitch")){			
			SxRNAEnvironment.isTranslationalSwitch=Boolean.parseBoolean(properties.getProperty("isTranslationalSwitch"));						
		}	
		if(properties.containsKey("foldCommand")){
			SxRNAEnvironment.foldCommand=properties.getProperty("foldCommand");						
		}
		if(properties.containsKey("cofoldCommand")){
			SxRNAEnvironment.cofoldCommand=properties.getProperty("cofoldCommand");						
		}		
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		StringBuffer sb = new StringBuffer();
		sb.append(super.toString());
		sb.append("SxRNAEnvironment [getFoldCommand()="+getFoldCommand());
		sb.append(", getCofoldCommand()="+getCofoldCommand());
		sb.append(", isTranslationalSwitch()="+isTranslationalSwitch);
		sb.append(", getTarget_trigger()="+getTarget_trigger());
		sb.append("]");
		return sb.toString();
	}

	
	
	


}
