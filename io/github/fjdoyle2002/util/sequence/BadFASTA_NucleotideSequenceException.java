package io.github.fjdoyle2002.util.sequence;

/**
 * Exception to denote occurence of an improper FASTA sequence
 */

public class BadFASTA_NucleotideSequenceException extends Exception
{
	public BadFASTA_NucleotideSequenceException(String gripe)
	{
		super(gripe);
	}
}
