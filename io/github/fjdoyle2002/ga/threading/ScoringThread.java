package io.github.fjdoyle2002.ga.threading;

import io.github.fjdoyle2002.ga.Environment;
import io.github.fjdoyle2002.ga.FitnessAlgorithm;
import io.github.fjdoyle2002.ga.Individual;

public class ScoringThread implements Runnable {
	
	private long start_time = 0;
	private Environment environment = null;
	private boolean should_run = true;
	private Individual current_individual = null;
	private Thread thread = null;
	private long id = -1;
	/**
	 * 
	 * @param fa
	 */
	public ScoringThread(Environment env){
		this.environment = env;
		this.thread = new Thread(this);
		this.id = this.thread.getId();
		this.thread.setDaemon(true);
		this.thread.start();		
	}
	/**
	 * 
	 */
	public void reset(){
		this.current_individual = null;
		this.start_time = System.currentTimeMillis();		
	}
	/**
	 * 
	 * @param should_run
	 */
	public void setShouldRun(boolean should_run){
		this.should_run  = should_run;
		
	}

	/**
	 * 
	 * @param current_individual
	 */
	public void scoreNewIndividual(Individual current_individual){
		this.current_individual = current_individual;
		this.start_time = System.currentTimeMillis();
	}
	/**
	 * 
	 * @return
	 */
	public boolean isReadyForNewIndividual(){
		return (this.should_run&&(this.current_individual==null || !this.current_individual.needsScoring()));
	}
	
	/**
	 * 
	 * @return
	 */
	public long getStartTime(){
		return this.start_time;
	}
	/**
	 * 
	 * @return
	 */
	public Individual getCurrentIndividual(){
		return this.current_individual;
	}
	

	
	
	@Override
	public void run() {
		while(this.should_run){
			try{
				Thread.sleep(100);
			}catch(InterruptedException ie){
				Thread.currentThread().interrupt();
				System.err.println("ScoringThread was interrupted with: "+this.current_individual);
			}
			if(this.current_individual != null && this.current_individual.needsScoring()){
				//update in case of delay 
				this.start_time=System.currentTimeMillis();
				this.environment.getFitness_algorithm().score(this.current_individual, this.id);
				long last_score_time = System.currentTimeMillis()-this.start_time;
				this.environment.registerScoreTime(last_score_time);
				
			}
			
		}
		
	}
	
	
}
