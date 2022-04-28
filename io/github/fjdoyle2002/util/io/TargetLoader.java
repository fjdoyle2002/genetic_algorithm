package io.github.fjdoyle2002.util.io;

import io.github.fjdoyle2002.util.TargetRNA;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


public class TargetLoader {
	/**
	 * Reads target RNA's from file (file must be single line... tab separated id and double precision value, no comments
	 */
	public static List<TargetRNA> loadTarget_RNAs(File targetRNAFile) {
		System.out.println("Reading target RNAs from file...");
		ArrayList<TargetRNA> ncRNA_list = new ArrayList<TargetRNA>();
		try{
			BufferedReader targetRNAReader = new BufferedReader(new FileReader(targetRNAFile));	
			String lineIn = targetRNAReader.readLine();

			while(lineIn!=null){
				if(!lineIn.startsWith("#") && !lineIn.trim().equals("")){
					String[] tokens = lineIn.split("\\t");
					String currentId = tokens[0].trim();
					double curr_off_val = 1;
					double curr_on_val = 1;
					if(tokens.length>1){
						curr_off_val = Double.parseDouble(tokens[1]);						
					}
					if(tokens.length>2){
						curr_on_val = Double.parseDouble(tokens[2]);						
					}
								
					ncRNA_list.add(new TargetRNA(currentId, curr_off_val, curr_on_val));
				}
				lineIn = targetRNAReader.readLine();				
			}
			targetRNAReader.close();
			//System.out.println("Finished reading ncRNAs from file.");
			//System.out.println("Number of ncRNAs loaded: "+ncRNA_list.size());
		}catch(IOException ioe){
			ioe.printStackTrace();
			System.err.println("IO Exception in loadTarget_RNAs for "+targetRNAFile.getName());
			System.exit(1);//don't attempt recovery on IO exception
		}
		return ncRNA_list;
		
	}
	
	public static Map<String,TargetRNA> loadTarget_RNA_as_map(File targetRNAFile) {
		System.out.println("Reading target RNAs from file...");
		HashMap<String, TargetRNA> ncRNA_map = new HashMap<String, TargetRNA>();
		try{
			BufferedReader targetRNAReader = new BufferedReader(new FileReader(targetRNAFile));	
			String lineIn = targetRNAReader.readLine();

			while(lineIn!=null){
				if(!lineIn.startsWith("#") && !lineIn.trim().equals("")){
					String[] tokens = lineIn.split("\\t");
					String currentId = tokens[0].trim();
					double curr_off_val = 1;
					double curr_on_val = 1;
					if(tokens.length>1){
						curr_off_val = Double.parseDouble(tokens[1]);						
					}
					if(tokens.length>2){
						curr_on_val = Double.parseDouble(tokens[2]);						
					}
					TargetRNA existing_entry = ncRNA_map.get(currentId);
					if(null != existing_entry){
						existing_entry.setOn_state_value(existing_entry.getOn_state_value()+curr_on_val);
						existing_entry.setOff_state_value(existing_entry.getOff_state_value()+curr_off_val);
						System.err.println("Found duplicate target data: "+currentId);
					}else{
						ncRNA_map.put(currentId, new TargetRNA(currentId, curr_off_val, curr_on_val));
					}
				}
				lineIn = targetRNAReader.readLine();				
			}
			targetRNAReader.close();
			//System.out.println("Finished reading ncRNAs from file.");
			//System.out.println("Number of ncRNAs loaded: "+ncRNA_list.size());
		}catch(IOException ioe){
			ioe.printStackTrace();
			System.err.println("IO Exception in loadTarget_RNAs for "+targetRNAFile.getName());
			System.exit(1);//don't attempt recovery on IO exception
		}
		return ncRNA_map;
		
	}

}
