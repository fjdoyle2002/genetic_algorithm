package io.github.fjdoyle2002.ga;

public abstract class Gene implements Cloneable{

	/**
	 * Constructor
	 */
	public Gene() {	}
	/**
	 * Mutation method
	 */
	public abstract void mutate();
	/* (non-Javadoc)
	 * @see java.lang.Object#clone()
	 */
	@Override
	protected Object clone() throws CloneNotSupportedException {
		// TODO Auto-generated method stub
		return super.clone();
	}
	
	

}
