package io.github.fjdoyle2002.util.sequence;


/**
 * @author Frank Doyle
 * @version 1.0
 */

/** A class to handle generic sequence manipulation and analysis requests.
 */

import java.util.*;

public class SequenceTools{
	protected static final String validBases = "[AaCcGgTtUu]+";//regular expression to match sequences of only specific nucleotide base codes
	protected static final String validFASTA = "[AaCcGgTtUuRrYyKkMmSsWwBbDdHhVvNn-]+";//regular expression to match sequences of only nucleotide/nucleotide-variable FASTA codes

//------------------------------------------------------------------------------------------------------
	public static void main(String[] args){
		//for testing only
		SequenceTools myToolSet = new SequenceTools();
		if(args.length < 1){
			System.out.println("Please provide arguments:\n    Space separated nucelotide base/FASTA variable - sequences.");
		}

		else{
			int i = 0;
			while (i < args.length){
				if(isValidNucleotideSequence(args[i])){
					try{
						System.out.println("gc content of:\n "+args[i]+"\nis: ");
						System.out.println("    "+myToolSet.gcContentOf(args[i]));
						System.out.println("Reverse complement of:\n "+args[i]+"\nis: ");
						System.out.println("    "+myToolSet.reverseComplementOf(args[i]));
						double[] composition = myToolSet.nucleotideCompostionOf(args[i]);
						System.out.println("Nucleotide percent composition (as {a,c,g,t(u)]:\n["+composition[0]+","+composition[1]+","+composition[2]+","+composition[3]+"]");

					}
					catch(Exception e){
						System.out.println("Error:\n  "+e);
					}
				}//end if(isValidNuc....
				else if(isValidFASTA_NucleotideSequence(args[i])){
					try{
						System.out.println("Minimum possible gc content of:\n "+ args[i]+"\nis: ");
						System.out.println("    "+myToolSet.min_gcContentInFasta(args[i]));
						System.out.println("Maximum possible gc content of:\n "+ args[i]+"\nis: ");
						System.out.println("    "+myToolSet.max_gcContentInFasta(args[i]));
						System.out.println("Minimum possible at content of:\n "+ args[i]+"\nis: ");
						System.out.println("    "+myToolSet.min_atContentInFasta(args[i]));
						System.out.println("Maximum possible at content of:\n "+ args[i]+"\nis: ");
						System.out.println("    "+myToolSet.max_atContentInFasta(args[i]));
						System.out.println("Reverse complement of:\n "+ args[i]+"\n is: ");
						System.out.println("    "+myToolSet.reverseComplementOf(args[i]));
					}
					catch(Exception e){
						System.out.println("Error:\n  "+e);
					}
				}//end if(isValidFAS....
				else{
					System.out.println(" Not a valid nucleotide base or nucleotide FASTA sequence: " + args[i]);
				}
				i++;
			}//end while
		}//end else
	}
//------------------------------------------------------------------------------------------------------
/*
  Used to calculate GC percentage in a given nucleotide string.  Accepts all valid FASTA nucleotide codes.
  However, only those with definitive GC or AT relationships are used as part of calculation.
  For instance sequence "agn"'s gc content would be 1/2 or 50%
  For max and minimum possible gc content of a FASTA nucleotide sequence, use min_gcContentInFasta() & max_gcContentInFasta()
  @param seq a nucleotide sequence to calculate GC content for.
  @return content The GC percentage of the given sequence
*/

	public static float gcContentOf(String seq) throws BadFASTA_NucleotideSequenceException{

		float content = 0;
		if(null==seq ||(!(seq.matches(validFASTA)))){
			throw new BadFASTA_NucleotideSequenceException("Improper sequence given: \""+ seq +"\".  Only valid FASTA nucleotide codes allowed (excluding - 'gap').");
		}//end if
		else if(seq.length() > 0){
			String bigCaseSequence = seq.toUpperCase();
			//System.out.println("Big case sequence is: "+bigCaseSequence);
			float gcCount = 0;//total g or c found in seq
			int i = 0;//position counter
			int baseCount = 0;
			while(i < bigCaseSequence.length()){
				switch (bigCaseSequence.charAt(i)){
					//only increase GC count for codes with definitive GC correlation
					//only increase base count for codes with definitive AT or GC correlation
					case 'C':
						gcCount++;
						baseCount++;
						break;
					case 'G':
						gcCount++;
						baseCount++;
						break;
					case 'S':
						gcCount++;
						baseCount++;
						break;
					case 'A':
						baseCount++;
						break;
					case 'T':
						baseCount++;
						break;
					case 'U':
						baseCount++;
						break;
					case 'W':
						baseCount++;
						break;
				}
				i++;
			}//end while
			//System.out.println("gcCount is: "+gcCount);
			//System.out.println("Big case sequence length is: "+bigCaseSequence.length())
			if(baseCount>0){
				content = gcCount/((float)baseCount)*100;
			}
			else{
				content = (float)0;
			}
		}//end else if
		return content;//content may be 0 at this stage, which is valid for 0 length sequence

	}//end gcContentOf(
	//------------------------------------------------------------------------------------------------------
	/*
	  Used to calculate nucleotide composition as percentages, returned in double[] as [a,c,g,t(u)]
	*/

		public static double[] nucleotideCompostionOf(String seq) throws BadFASTA_NucleotideSequenceException{

			int[] compositionCount = new int[4];
			double[] compositionPercent = new double[4];
			if(null==seq ||(!(seq.matches(validFASTA)))){
				throw new BadFASTA_NucleotideSequenceException("Improper sequence given: \""+ seq +"\".  Only valid FASTA nucleotide codes allowed (excluding - 'gap').");
			}//end if
			else if(seq.length() > 0){
				String bigCaseSequence = seq.toUpperCase();
				//System.out.println("Big case sequence is: "+bigCaseSequence);
				
				int i = 0;//position counter
				int baseCount = 0;
				while(i < bigCaseSequence.length()){
					switch (bigCaseSequence.charAt(i)){
						//only increase GC count for codes with definitive GC correlation
						//only increase base count for codes with definitive AT or GC correlation
						case 'C':
							compositionCount[1]++;
							baseCount++;
							break;
						case 'G':
							compositionCount[2]++;
							baseCount++;
							break;
						case 'S':
							baseCount++;
							break;
						case 'A':
							compositionCount[0]++;
							baseCount++;
							break;
						case 'T':
							compositionCount[3]++;
							baseCount++;
							break;
						case 'U':
							compositionCount[3]++;
							baseCount++;
							break;
						case 'W':
							baseCount++;
							break;
					}
					i++;
				}//end while
				//System.out.println("gcCount is: "+gcCount);
				//System.out.println("Big case sequence length is: "+bigCaseSequence.length())
				if(baseCount>0){
					for(int pos=0;pos<4;pos++){
						compositionPercent[pos] = (float)compositionCount[pos]/((float)baseCount)*100d;
					}
				}

			}//end else if
			return compositionPercent;//content may be 0 at this stage, which is valid for 0 length sequence

		}//end gcContentOf(
//------------------------------------------------------------------------------------------------------
/*
  Used to calculate percentage of N or n characters in a given nucleotide string.  Accepts all valid FASTA nucleotide codes.
  @param seq a nucleotide sequence to calculate N content for.
  @return content The N|n percentage of the given sequence
*/

	public static float nContentOf(String seq) throws BadFASTA_NucleotideSequenceException{

		float content = 0;
		if(null==seq ||(!(seq.matches(validFASTA)))){
			throw new BadFASTA_NucleotideSequenceException("Improper sequence given: \""+ seq +"\".  Only valid FASTA nucleotide codes allowed (excluding - 'gap').");
		}//end if
		else if(seq.length() > 0){
			String bigCaseSequence = seq.toUpperCase();
			//System.out.println("Big case sequence is: "+bigCaseSequence);
			float nCount = 0;//total g or c found in seq
			int i = 0;//position counter
			int baseCount = 0;
			while(i < bigCaseSequence.length()){
				switch (bigCaseSequence.charAt(i)){
					//only increase GC count for codes with definitive GC correlation
					//only increase base count for codes with definitive AT or GC correlation
					case 'N':
						nCount++;
						baseCount++;
						break;
					case 'A':
						baseCount++;
						break;
					case 'C':
						baseCount++;
						break;
					case 'G':
						baseCount++;
						break;
					case 'T':
						baseCount++;
						break;
					case 'U':
						baseCount++;
						break;
					case 'R':
						baseCount++;
						break;
					case 'Y':
						baseCount++;
						break;
					case 'K':
						baseCount++;
						break;
					case 'M':
						baseCount++;
						break;
					case 'S':
						baseCount++;
						break;
					case 'W':
						baseCount++;
						break;
					case 'B':
						baseCount++;
						break;
					case 'D':
						baseCount++;
						break;
					case 'H':
						baseCount++;
						break;
					case 'V':
						baseCount++;
						break;

				}
				i++;
			}//end while
			//System.out.println("gcCount is: "+gcCount);
			//System.out.println("Big case sequence length is: "+bigCaseSequence.length())
			if(baseCount>0){
				content = nCount/((float)baseCount)*100;
			}
			else{
				content = (float)0;
			}
		}//end else if
		return content;//content may be 0 at this stage, which is valid for 0 length sequence or sequence with no GC content

	}//end gcContentOf(
//------------------------------------------------------------------------------------------------------
/*
  Used to calculate AT percentage in a given nucleotide string.  Accepts both U and T as valid bases.
  For max and minimum possible AT content of a FASTA nucleotide sequence, use min_atContentInFasta() & max_atContentInFasta()
  @param seq a nucleotide sequence to calculate AT content for.
  @return content The AT percentage of the given sequence
*/

	public static final float atContentOf(String seq) throws BadFASTA_NucleotideSequenceException{
		return (float)100 - gcContentOf(seq);
	}
//------------------------------------------------------------------------------------------------------
/*
  Used to calculate minimum possible GC percentage in a given FASTA nucleotide string.
  @param seq a FASTA nucleotide sequence to calculate minimum GC content for.
  @return content The minimum GC percentage of the given sequence
*/

	public static float min_gcContentInFasta(String seq) throws BadFASTA_NucleotideSequenceException{

		float minContent = 0;
		if(null==seq ||(!(seq.matches(validFASTA)))){
			throw new BadFASTA_NucleotideSequenceException("Improper fasta nucleotide sequence given: \""+ seq+"\".");
		}//end if
		else if(seq.length() > 0){
			String bigCaseSequence = seq.toUpperCase();
			//System.out.println("Big case sequence is: "+bigCaseSequence);
			float gcCount = 0;//total g or c found in seq
			int i = 0;//position counter
			int baseCount = bigCaseSequence.length();
			while(i < bigCaseSequence.length()){
				switch (bigCaseSequence.charAt(i)){
					case 'C':
						gcCount++;
						break;
					case 'c':
						gcCount++;
						break;
					case 'G':
						gcCount++;
					case 'g':
						gcCount++;
						break;
					case 'S':
						gcCount++;
						break;
					case 's':
						gcCount++;
						break;
					case '-':
						baseCount--;//do not use gap positions in calculation
						break;

				}//end switch
				i++;
			}//end while
			//System.out.println("gcCount is: "+gcCount);
			//System.out.println("Big case sequence length is: "+bigCaseSequence.length());
			if (baseCount > 0){
				minContent = gcCount/((float)baseCount)*100;
			}
			else{
				minContent = (float)0;
			}
		}//end else if
		return minContent;//content may be 0 at this stage, which is valid for 0 length sequence

	}//end min_gcContentInFasta(
//------------------------------------------------------------------------------------------------------
/*
  Used to calculate maximum possible GC percentage in a given FASTA nucleotide string.
  @param seq a FASTA nucleotide sequence to calculate maximum GC content for.
  @return content The maximum GC percentage of the given sequence
*/

	public static float max_gcContentInFasta(String seq) throws BadFASTA_NucleotideSequenceException{

		float maxContent = 0;
		if(null==seq ||(!(seq.matches(validFASTA)))){
			throw new BadFASTA_NucleotideSequenceException("Improper fasta nucleotide sequence given: "+ seq);
		}//end if
		else if(seq.length() > 0){
			String bigCaseSequence = seq.toUpperCase();
			//System.out.println("Big case sequence is: "+bigCaseSequence);
			float gcCount = 0;//total g or c found in seq
			int i = 0;//position counter
			while(i < bigCaseSequence.length()){
				switch (bigCaseSequence.charAt(i)){
					case 'C':
						gcCount++;
						break;
					case 'G':
						gcCount++;
						break;
					case 'R':
						gcCount++;
						break;
					case 'Y':
						gcCount++;
						break;
					case 'K':
						gcCount++;
						break;
					case 'M':
						gcCount++;
						break;
					case 'S':
						gcCount++;
						break;
					case 'B':
						gcCount++;
						break;
					case 'D':
						gcCount++;
						break;
					case 'H':
						gcCount++;
						break;
					case 'V':
						gcCount++;
						break;
					case 'N':
						gcCount++;
						break;
				}//end switch
				i++;
			}//end while
			//System.out.println("gcCount is: "+gcCount);
			//System.out.println("Big case sequence length is: "+bigCaseSequence.length());
			maxContent = gcCount/((float)seq.length())*100;
		}//end else if
		return maxContent;//content may be 0 at this stage, which is valid for 0 length sequence

	}//end max_gcContentInFasta(

//------------------------------------------------------------------------------------------------------
/*
  Used to calculate minimum possible AT percentage in a given FASTA nucleotide string.
  @param seq a FASTA nucleotide sequence to calculate minimum AT content for.
  @return content The minimum AT percentage of the given sequence
*/

	public static final float min_atContentInFasta(String seq) throws BadFASTA_NucleotideSequenceException{
		return (float)100 - max_gcContentInFasta(seq);
	}
//------------------------------------------------------------------------------------------------------
/*
  Used to calculate maximum possible AT percentage in a given FASTA nucleotide string.
  @param seq a FASTA nucleotide sequence to calculate maximum AT content for.
  @return content The maximum AT percentage of the given sequence
*/

	public static final float max_atContentInFasta(String seq) throws BadFASTA_NucleotideSequenceException{
		return (float)100 - min_gcContentInFasta(seq);
	}
//------------------------------------------------------------------------------------------------------
/*
  Used to get the complement of a given nucleotide sequence(gcta returns cgat).
  @param seq a nucleotide sequence to create complement of
  @return newSeq The complement of the given sequence
*/

	public static String complementOf(String seq) throws BadFASTA_NucleotideSequenceException{

		String newSeq = "";

		if(null==seq ||(!(seq.matches(validFASTA)))){
			throw new BadFASTA_NucleotideSequenceException("Improper nucleotide sequence given: "+ seq +" .Allowable bases: AaCcGgTtUuNn");
		}//end if
		else if(seq.length() > 0){
			String bigCaseSequence = seq.toUpperCase();
			//System.out.println("Big case sequence is: "+bigCaseSequence);
			StringBuffer newSequenceBuffer = new StringBuffer();
			int i = 0;//position counter
			while(i < bigCaseSequence.length()){
				switch (bigCaseSequence.charAt(i)){
					case 'A':
						newSequenceBuffer.append("T");
						break;
					case 'C':
						newSequenceBuffer.append("G");
						break;
					case 'G':
						newSequenceBuffer.append("C");
						break;
					case 'T':
						newSequenceBuffer.append("A");
						break;
					case 'U':
						newSequenceBuffer.append("A");
						break;
					case 'R':
						newSequenceBuffer.append("Y");
						break;
					case 'Y':
						newSequenceBuffer.append("R");
						break;
					case 'K':
						newSequenceBuffer.append("M");
						break;
					case 'M':
						newSequenceBuffer.append("K");
						break;
					case 'S':
						newSequenceBuffer.append("W");
						break;
					case 'W':
						newSequenceBuffer.append("S");
						break;
					case 'B':
						newSequenceBuffer.append("V");
						break;
					case 'D':
						newSequenceBuffer.append("H");
						break;
					case 'H':
						newSequenceBuffer.append("D");
						break;
					case 'V':
						newSequenceBuffer.append("B");
						break;
					case 'N':
						newSequenceBuffer.append("N");
						break;
				}
				i++;
			}//end while
			//System.out.println("gcCount is: "+gcCount);
			//System.out.println("Big case sequence length is: "+bigCaseSequence.length());
			newSeq = newSequenceBuffer.toString();
		}//end else if
		return newSeq;

	}//end complementOf(
//------------------------------------------------------------------------------------------------------
/*
  Used to get the reverse complement of a given nucleotide sequence(gcta returns tagc).
  @param seq a nucleotide sequence to create reverse complement of
  @return newSeq The reverse complement of the given sequence
*/

	public static String reverseComplementOf(String seq) throws BadFASTA_NucleotideSequenceException{
		return reverseComplementOf(seq, false);
	}

//------------------------------------------------------------------------------------------------------
/*
  Used to get the reverse complement of a given nucleotide sequence(gcta returns tagc).
  @param seq a nucleotide sequence to create reverse complement of
  @return newSeq The reverse complement of the given sequence
*/

	public static String reverseComplementOf(String seq, boolean convertThymidine) throws BadFASTA_NucleotideSequenceException{

		String newSeq = "";

		if(null==seq ||(!(seq.matches(validFASTA)))){
			throw new BadFASTA_NucleotideSequenceException("Improper nucleotide sequence given: "+ seq +" .Allowable bases: AaCcGgTtUuNn");
		}//end if
		else if(seq.length() > 0){
			//String bigCaseSequence = seq.toUpperCase();
			//System.out.println("Big case sequence is: "+bigCaseSequence);
			StringBuffer newSequenceBuffer = new StringBuffer();
			int i = (seq.length()-1);//position counter
			while(i >= 0){
				switch (seq.charAt(i)){
					case 'A':
						if(convertThymidine)
							newSequenceBuffer.append("U");
						else
							newSequenceBuffer.append("T");
						break;
					case 'a':
						if(convertThymidine)
							newSequenceBuffer.append("u");
						else
							newSequenceBuffer.append("t");
						break;
					case 'C':
						newSequenceBuffer.append("G");
						break;
					case 'c':
						newSequenceBuffer.append("g");
						break;
					case 'G':
						newSequenceBuffer.append("C");
						break;
					case 'g':
						newSequenceBuffer.append("c");
						break;
					case 'T':
						newSequenceBuffer.append("A");
						break;
					case 't':
						newSequenceBuffer.append("a");
						break;
					case 'U':
						newSequenceBuffer.append("A");
						break;
					case 'u':
						newSequenceBuffer.append("a");
						break;
					case 'R':
						newSequenceBuffer.append("Y");
						break;
					case 'r':
						newSequenceBuffer.append("y");
						break;
					case 'Y':
						newSequenceBuffer.append("R");
						break;
					case 'y':
						newSequenceBuffer.append("r");
						break;
					case 'K':
						newSequenceBuffer.append("M");
						break;
					case 'k':
						newSequenceBuffer.append("m");
						break;
					case 'M':
						newSequenceBuffer.append("K");
						break;
					case 'm':
						newSequenceBuffer.append("k");
						break;
					case 'S':
						newSequenceBuffer.append("W");
						break;
					case 's':
						newSequenceBuffer.append("w");
						break;
					case 'W':
						newSequenceBuffer.append("S");
						break;
					case 'w':
						newSequenceBuffer.append("s");
						break;
					case 'B':
						newSequenceBuffer.append("V");
						break;
					case 'b':
						newSequenceBuffer.append("v");
						break;
					case 'D':
						newSequenceBuffer.append("H");
						break;
					case 'd':
						newSequenceBuffer.append("h");
						break;
					case 'H':
						newSequenceBuffer.append("D");
						break;
					case 'h':
						newSequenceBuffer.append("d");
						break;
					case 'V':
						newSequenceBuffer.append("B");
						break;
					case 'v':
						newSequenceBuffer.append("b");
						break;
					case 'N':
						newSequenceBuffer.append("N");
						break;
					case 'n':
						newSequenceBuffer.append("n");
						break;
					case '-':
						newSequenceBuffer.append("-");
						break;
				}
				i--;
			}//end while
			//System.out.println("gcCount is: "+gcCount);
			//System.out.println("Big case sequence length is: "+bigCaseSequence.length());
			newSeq = newSequenceBuffer.toString();
		}//end else if
		return newSeq;

	}//end reverseComplementOf(
	//------------------------------------------------------------------------------------------------------
	/*
	  Used to get the reverse of a given nucleotide sequence(gcta returns atcg).
	  @param seq a nucleotide sequence to create reverse complement of
	  @return newSeq The reverse of the given sequence
	*/

		public static String reverseOf(String seq) throws BadFASTA_NucleotideSequenceException{

			String newSeq = "";

			if(null==seq ||(!(seq.matches(validFASTA)))){
				throw new BadFASTA_NucleotideSequenceException("Improper nucleotide sequence given: "+ seq +" .Allowable bases: AaCcGgTtUuNn");
			}//end if
			else if(seq.length() > 0){
				//String bigCaseSequence = seq.toUpperCase();
				//System.out.println("Big case sequence is: "+bigCaseSequence);
				StringBuffer newSequenceBuffer = new StringBuffer();
				int i = (seq.length()-1);//position counter
				while(i >= 0){
					newSequenceBuffer.append(seq.charAt(i));
					i--;
				}
				newSeq = newSequenceBuffer.toString();
			}//end else if
			return newSeq;

		}//end reverseOf(	
//------------------------------------------------------------------------------------------------------
/*
  Used to check a sequence for the appearance of repeat masked (lower case) regions.
  @param seq a FASTA nucleotide sequence to check for repeats
  @return hasRepeats The boolean value telling whether this sequence contains lowercase, repeat designators
*/

	public static final boolean repeatsMarkedIn(String seq){
		boolean hasRepeats = false;
		if(!(null == seq)){
			String bigCaseSequence = seq.toUpperCase();
			hasRepeats = (!(seq.equals(bigCaseSequence)));
		}
		return hasRepeats;
	}
//------------------------------------------------------------------------------------------------------
/*
  Used to test for proper nucleotide sequence format (ie, contains only acgt or u, case insensitive).
  0 length sequences are considered to be invalid.
  @param seq a nucleotide sequence to test
  @return isValid True or false
*/

	public static boolean isValidNucleotideSequence(String seq){
		boolean isValid;
		if(null==seq ||(seq.length() == 0)){
			isValid = false;
		}//end if
		else if(seq.matches(validBases)){
			isValid = true;
		}
		else{
			isValid = false;
		}
		return isValid;
	}
/*
  Used to test for proper FASTA nucleotide sequence format (ie, contains only FASTA standard nucleotide codes, not amino acids ).
  0 length sequences are considered to be invalid.
  @param seq a nucleotide sequence to test
  @return isValid True or false
*/

	public static boolean isValidFASTA_NucleotideSequence(String seq){
		boolean isValid;
		if(null==seq ||(seq.length() == 0)){
			isValid = false;
		}//end if
		else if(seq.matches(validFASTA)){
			isValid = true;
		}
		else{
			isValid = false;
		}
		return isValid;
	}

}//end class sequenceTools

