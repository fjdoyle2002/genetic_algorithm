package io.github.fjdoyle2002.parsers;
import java.util.Stack;
/**
 * @author fd299212
 * Class provides logic for parsing "dot parentheses" designation of 
 * secondary structure as produced by RNAcofold  
 * 
 */
public class DotParenthesesParser {


	public static int[] parse(String dotParenthesesLine){
		//System.out.println("DotParenthesesParser.parse() entered with:\n\t"+dotParenthesesLine);
		Stack<Integer> leftParenthesesPositions = null;
		int[] matchPositions = null;
		if(!(null==dotParenthesesLine) && dotParenthesesLine.length()>0)
		leftParenthesesPositions = new Stack<Integer>();
		matchPositions=new int[dotParenthesesLine.length()];
		for(int i=0;i<dotParenthesesLine.length();i++){
			char c=dotParenthesesLine.charAt(i);
			switch (c){
				case '(':
					leftParenthesesPositions.push(i);
					break;
				case '.':
					matchPositions[i]=-1;
					break;
				case ')':
					int leftMatchPos=leftParenthesesPositions.pop();
					matchPositions[leftMatchPos]=i;
					matchPositions[i]=leftMatchPos;
					break;
				case '&':
					matchPositions[i]=-2;
					break;			
				default:
					//throw exception
					System.err.println("Invalid Character for DotParenthesesParser: '"+c+"' in |"+dotParenthesesLine+"| at position:"+i  );				
			}
			
		}
		/*for(int i=0;i<matchPositions.length;i++){
			System.out.print(","+matchPositions[i]);
		}
		System.out.println();*/
		
		return matchPositions;
	}
	
	
}
