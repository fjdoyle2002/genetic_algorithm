package io.github.fjdoyle2002.util.io;
import io.github.fjdoyle2002.ga.Individual;
import io.github.fjdoyle2002.util.NonCodingRNA;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.Iterator;
import java.util.List;



public class SimpleLog {
		
		private PrintWriter logfile = null;
		private String filename = null;
		/**
		 * Constructor 
		 * @param logname
		 */
		public SimpleLog(String logname){
			filename = System.getProperty("user.dir")
					+ System.getProperty("file.separator") + logname;
			try {
				logfile = new PrintWriter(new FileWriter(filename, false), true);
			} catch (Exception e) {
				System.out.println("Exception in creating file: " +filename +"\n"+ e.getLocalizedMessage());
			}
		}
		/**
		 * Writes a line of text to the PrintWriter
		 * @param comment
		 */
		public synchronized void comment(String comment){
			if(logfile != null){
				logfile.println(comment);	
			}else{
				System.out.println(comment);
			}
		}
		/**
		 * 
		 * @param objList
		 */
		public synchronized void comment(List<NonCodingRNA> ncList){
			StringBuffer sb = new StringBuffer();
			sb.append("List:");
			if (ncList != null){
				Iterator<NonCodingRNA> ncIter = ncList.iterator();
				while(ncIter.hasNext()){
					sb.append(ncIter.next().toString());
					if(ncIter.hasNext()){
						sb.append(",");
					}
				}
				
			}
			if(this.logfile != null){
				comment(sb);
			}else{
				System.out.println(sb.toString());
			}
			
		}
		/**
		 * Writes a line of text to the PrintWriter
		 * @param comment
		 */
		public synchronized void comment(StringBuffer comment){
			if(logfile != null){
				logfile.println(comment.toString());	
			}else{
				System.out.println(comment.toString());
			}
		}
		/**
		 * Closes the PrintWriter
		 */
		public synchronized void closeLogfile(){
			try{
				logfile.flush();
				logfile.close();
			}catch(Exception e){
				System.out.println("Exception closing logfile: "+filename+"\n"+e.getLocalizedMessage());
			}
		}
		/**
		 * 
		 * @param individual
		 */
		public synchronized void comment(Individual individual) {
			if(logfile != null){
				logfile.println(individual.toString());	
			}else{
				System.out.println(individual.toString());
			}
			
		}

}
