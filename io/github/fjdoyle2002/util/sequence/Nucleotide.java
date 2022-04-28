package io.github.fjdoyle2002.util.sequence;

public interface Nucleotide {
	public static final int A = 0;
	public static final int C = 1;
	public static final int G = 2;
	public static final int T = 3;
	public static final int U = 4;
	
	public static final char[] bases = {'A','C','G','T','U'};
	
	public static final int[] Y = {C,T,U};
	public static final int[] R = {G,A};
	public static final int[] H = {A,C,T,U};
	public static final int[] K = {G,T,U};
	public static final int[] S = {G,C};
	public static final int[] M = {A,C};
	public static final int[] W = {A,T,U};
	public static final int[] B = {G,T,U,C};
	public static final int[] D = {G,A,T,U};
	public static final int[] V = {G,C,A};
	public static final int[] N = {A,C,G,T,U};
	
	
	
}
