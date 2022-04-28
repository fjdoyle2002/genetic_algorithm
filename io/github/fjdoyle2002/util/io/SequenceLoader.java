package io.github.fjdoyle2002.util.io;

import io.github.fjdoyle2002.util.NonCodingRNA;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;




public class SequenceLoader {


	/**
	 * Reads nc RNA's from file (file must be single line header, followed by single line seq, no space, no comments
	 */
	public static List<NonCodingRNA> loadNC_RNAs_as_list(File ncRNAFile) {
		//System.out.println("Reading ncRNAs from file...");
		ArrayList<NonCodingRNA> ncRNA_list = new ArrayList<NonCodingRNA>();
		try{
			BufferedReader ncRNAReader = new BufferedReader(new FileReader(ncRNAFile));	
			String lineIn = ncRNAReader.readLine();
			String currentHeader = null;
			String currentSeq = null;

			while(lineIn!=null){
				if(lineIn.startsWith(">")){
					currentHeader = lineIn.replace(">", "");
				}else if(!(lineIn.trim().equals("")||lineIn.startsWith("#"))){
					currentSeq = lineIn.trim();
					if(!(null==currentHeader||null==currentSeq)){
						ncRNA_list.add(new NonCodingRNA(currentHeader, currentSeq));
						//full_ncRNA_map.put(currentHeader, new NonCodingRNA(currentHeader, currentSeq));
					}
					currentHeader = null;
					currentSeq = null;
				}
				lineIn = ncRNAReader.readLine();
				
			}
			//System.out.println("Finished reading ncRNAs from file.");
			//System.out.println("Number of ncRNAs loaded: "+ncRNA_list.size());
		}catch(IOException ioe){
			ioe.printStackTrace();
			System.err.println("IO Exception in loadNC_RNAs for "+ncRNAFile.getName());
			System.exit(1);//don't attempt recovery on IO exception
		}
		return ncRNA_list;
		
	}
	
	/**
	 * Reads nc RNA's from file (file must be single line header, followed by single line seq, no space, no comments
	 */
	public static Map<String,NonCodingRNA> loadNC_RNAs_as_map(File ncRNAFile) {
		//System.out.println("Reading ncRNAs from file...");
		HashMap<String, NonCodingRNA> ncRNA_map = new HashMap<String, NonCodingRNA>();
		try{
			BufferedReader ncRNAReader = new BufferedReader(new FileReader(ncRNAFile));	
			String lineIn = ncRNAReader.readLine();
			String currentHeader = null;
			String currentSeq = null;

			while(lineIn!=null){
				if(lineIn.startsWith(">")){
					currentHeader = lineIn.replace(">", "");
				}else if(!(lineIn.trim().equals("")||lineIn.startsWith("#"))){
					currentSeq = lineIn.trim();
					if(!(null==currentHeader||null==currentSeq)){
						NonCodingRNA currNcRNA = new NonCodingRNA(currentHeader, currentSeq); 
						ncRNA_map.put(currNcRNA.getId(), currNcRNA);
						//full_ncRNA_map.put(currentHeader, new NonCodingRNA(currentHeader, currentSeq));
					}
					currentHeader = null;
					currentSeq = null;
				}
				lineIn = ncRNAReader.readLine();
				
			}
			//System.out.println("Finished reading ncRNAs from file.");
			//System.out.println("Number of ncRNAs loaded: "+ncRNA_list.size());
		}catch(IOException ioe){
			ioe.printStackTrace();
			System.err.println("IO Exception in loadNC_RNAs for "+ncRNAFile.getName());
			System.exit(1);//don't attempt recovery on IO exception
		}
		return ncRNA_map;
		
	}
}
