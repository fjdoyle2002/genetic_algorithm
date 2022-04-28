package io.github.fjdoyle2002.ga.threading;

import io.github.fjdoyle2002.ga.Environment;
import io.github.fjdoyle2002.ga.Individual;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;


public class ScoringThreadPool {
	private Environment environment = null;
	int maximum_threads = 1;
	List<ScoringThread> pool = null;

	
	public ScoringThreadPool(Environment environment){
		this.environment = environment;
		Environment.getLog().comment("Initializing ScoringThreadPool");		
		this.maximum_threads = this.environment.getMaximum_threads();
		int available_processors = Runtime.getRuntime().availableProcessors();
		//don't create more threads than available processors to use them
		if(this.maximum_threads > available_processors){
			Environment.getLog().comment("Adjusted requested thread number for available processors ["+this.maximum_threads+" to "+available_processors+"]");
			this.maximum_threads = available_processors;			
		}
		this.pool = new ArrayList<ScoringThread>();
		fillPool();
	}

	private void fillPool() {
		if(this.pool.size() < maximum_threads){	
			Environment.getLog().comment("ScoringThreadPool currently has ["+pool.size()+"] scoring threads");
			while(pool.size() < maximum_threads){
				pool.add(new ScoringThread(this.environment));
			}
			Environment.getLog().comment("Filled ScoringThreadPool with ["+pool.size()+"] scoring threads");
		}
	}
	/**
	 * 
	 * @param individualsToBeScored
	 */
	public void scoreIndividuals(List<Individual> individualsToBeScored){
		
		while(individualsToBeScored.size()>0){
			Iterator<Individual> individual_iter = individualsToBeScored.iterator();
			while(individual_iter.hasNext()){
				Individual current_individual = individual_iter.next();
				Iterator<ScoringThread> thread_iter = pool.iterator();
				while(thread_iter.hasNext()){
					ScoringThread currentThread = thread_iter.next();
					if(currentThread.isReadyForNewIndividual()){
						currentThread.scoreNewIndividual(current_individual);
						individual_iter.remove();
						break;
					}else{
						long runtime = System.currentTimeMillis() - currentThread.getStartTime();
						if(runtime > this.environment.getMaximum_time_to_wait_for_thread()){
							currentThread.setShouldRun(false);
							thread_iter.remove();
						}
					}
					
				}				
			}
			try{
				Thread.sleep(10);
			}catch(InterruptedException ie){
				Thread.currentThread().interrupt();
				System.err.println("rankMultiThread was interrupted with: "+ie.getLocalizedMessage());
			}
			//ensure we replenish scoring threads if any have been removed from pool
			//due to fatal conditions
			fillPool();
			
			
		}
		
	}
	/**
	 * 
	 */
	public void shutDown(){
		Iterator<ScoringThread> thread_iter = pool.iterator();
		while(thread_iter.hasNext()){
			ScoringThread currentThread = thread_iter.next();
			currentThread.setShouldRun(false);			
		}		
	}
		
}
