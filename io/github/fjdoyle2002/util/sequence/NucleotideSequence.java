package io.github.fjdoyle2002.util.sequence;

public abstract class NucleotideSequence {
	private StringBuffer sequence;
	private int type;
	public static final int GENOMIC = 0;
	public static final int MRNA = 1;
	public static final int _5PUTR = 2;
	public static final int ORF = 3;
	public static final int _3PUTR = 4;
	public static final int OLIGO =5;
	/**
	 * 
	 */
	public NucleotideSequence() {
		// TODO Auto-generated constructor stub
	}
	/**
	 * @return the sequence
	 */
	public StringBuffer getSequence() {
		return sequence;
	}
	/**
	 * @return the type
	 */
	public int getType() {
		return type;
	}
	

}
