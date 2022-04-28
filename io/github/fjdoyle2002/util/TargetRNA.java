package io.github.fjdoyle2002.util;

public class TargetRNA {
	private String id = null;
	private double off_state_value = 1;
	private double on_state_value = 1;


	
	
	public TargetRNA(String id, double off_state_value, double on_state_value) {
		super();
		this.id = id;
		this.off_state_value = off_state_value;
		this.on_state_value = on_state_value;

	}
	
	public String getId() {
		return id;
	}
	public void setId(String id) {
		this.id = id;
	}

	public double getOff_state_value() {
		return off_state_value;
	}

	public void setOff_state_value(double off_state_value) {
		this.off_state_value = off_state_value;
	}

	public double getOn_state_value() {
		return on_state_value;
	}

	public void setOn_state_value(double on_state_value) {
		this.on_state_value = on_state_value;
	}
	
	public double getPosTargetDifferentialValue(){
		return Math.abs(this.on_state_value - this.off_state_value);
	}
	
	public double getNegTargetDifferentialValue(){
		return Math.abs(this.off_state_value - this.on_state_value);
	}
	
	public double getAverageValue(){
		return (this.off_state_value + this.on_state_value)/2d;
	}
	
	public double getOnFoldChange(){
		double fold_change = 1;
		if(this.off_state_value>0){
			fold_change = this.on_state_value/this.off_state_value;
		}else{
			fold_change = this.on_state_value/.001;
		}
		return fold_change;
	}
	public double getOffFoldChange(){
		double fold_change = 1;
		if(this.on_state_value>0){
			fold_change = this.off_state_value/this.on_state_value;
		}else{
			fold_change = this.off_state_value/.001;
		}
		return fold_change;
	}	
	public double getAverageScore(){
		return ((this.off_state_value+this.on_state_value)/2d);
	} 


	@Override
	public String toString() {
		return "TargetRNA [id=" + id + ", off_state_value=" + off_state_value
				+ ", on_state_value=" + on_state_value + "]";
	}
	
	

}
