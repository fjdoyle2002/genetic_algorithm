package io.github.fjdoyle2002.util;

import io.github.fjdoyle2002.util.io.SequenceLoader;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.List;


public class NcRNAPicker {
	private String full_ncrna_filename = null;
	private String filter_ncrna_filename = null;
	private List<NonCodingRNA> ncrnas = null;
	 

	public NcRNAPicker(String[] args) {
		this.full_ncrna_filename = args[0];
		this.filter_ncrna_filename = args[1];
		
	}

	public static void main (String[] args){
		
		NcRNAPicker ncrnapicker = new NcRNAPicker(args);
		ncrnapicker.pick();
		
	}

	private void pick() {
		File full_ncrna_file = new File(this.full_ncrna_filename);
	    ncrnas = SequenceLoader.loadNC_RNAs_as_list(full_ncrna_file);
	    this.readFilterFile();

	}
	
	private void readFilterFile(){
		try{
			BufferedReader txReader = new BufferedReader(new FileReader(filter_ncrna_filename));	
			String lineIn = txReader.readLine();
			while(lineIn!=null){
				if(!(lineIn.trim().equals("")||lineIn.startsWith("#"))){

					String[] tokens = lineIn.split("\\t");
					String precursor = tokens[0].trim();
					String name = tokens[1].trim();
					Iterator<NonCodingRNA> nc_iter = this.ncrnas.iterator();
					while(nc_iter.hasNext()){
						NonCodingRNA curr = nc_iter.next();
						if(curr.getHeader().contains(name)){
							System.out.println(">"+curr.getHeader());
							System.out.println(curr.getSequence());
							break;
						}
					}
					
					
				}
				lineIn = txReader.readLine();				
			}
			txReader.close();
			System.out.println("Finished reading ncRNAs from file.");
			System.out.println("Number of ncRNAs loaded: "+ncrnas.size());
		}catch(IOException ioe){
			ioe.printStackTrace();
			System.err.println("IO Exception in loadNcRNA for "+filter_ncrna_filename);
			System.exit(1);//don't attempt recovery on IO exception
		}		
	}
	
	
}
