package io.github.fjdoyle2002.ga.sxrna;

import io.github.fjdoyle2002.ga.Gene;
import io.github.fjdoyle2002.ga.Individual;
import io.github.fjdoyle2002.util.NonCodingRNA;
import io.github.fjdoyle2002.util.Ribonucleotides;
import io.github.fjdoyle2002.util.SequenceTools;
import io.github.fjdoyle2002.util.TargetRNA;
import io.github.fjdoyle2002.util.sequence.BadFASTA_NucleotideSequenceException;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.ibm.compbio.seqalign.SmithWaterman;

public class SxRNAIndividual extends Individual {

	private static String sequence_nucleotide_constraints = null;
	private static String sequence_metadata_constraints = null;
	private static String motif_structure = null;
	private static List<String> motifs = null;

	private static int begin_absolute_motif_struct_pos = 0;
	private static int end_absolute_motif_struct_pos = 0;
	private int begin_motif_pos_in_fasta = 0;
	private int end_motif_pos_in_fasta = 0;
	private boolean trigger_spans_motif = false;
	private boolean motif_present_in_solo_fold = false;
	private boolean motif_present_in_cofold = false;
	
	private String rnaFold = null;
	private String rnaCofold = null;	
	
	private String nonTargetPosHits = "";
	private String nonTargetSuppressorHits = "";

	private String negTargetPosHits = "";
	private String negTargetSuppressorHits = "";
	
	private String posTargetPosHits = "";
	private String posTargetSuppressorHits = "";

	private double posTargetPosHitScore = 0;
	private double posTargetSuppressorHitScore = 0;
	private double nonTargetPosHitScore = 0;
	private double nonTargetSuppressorHitScore = 0;
	private double negTargetPosHitScore = 0;
	private double negTargetSuppressorHitScore = 0;

	/**
	 * Constructor
	 */
	public SxRNAIndividual(){
		super();

		createNucleotideGenes();

	}
	

	/**
	 * Sets the static variables needed to produce instances of this class.
	 * @param sequence_constraints
	 *        The FASTA code defining flexibility of nucleotide choice at each position. 
	 *        The length of this string determines the chromosome length for all sxRNAIndividuals. 
	 * @param sequence_metadata
	 *        Extensible code defining additional constraints beyond simple nucleotide choice.
	 *        For instance, whether or not a position in the chromosome can have a 'blank' allele.
	 * @param motif_structure
	 *        Dot parantheses notation defining the structure component of the desired motif
	 */
	public static void initializeStaticClassVariables(String sequence_nucleotide_constraints, String sequence_metadata_constrainst, String motif_structure)throws BadFASTA_NucleotideSequenceException{
		if(SequenceTools.isValidFASTA_NucleotideSequence(sequence_nucleotide_constraints)){
			SxRNAIndividual.sequence_nucleotide_constraints = sequence_nucleotide_constraints;
		}else{
			throw new BadFASTA_NucleotideSequenceException("sequence_nucleotide_constraints");
		}
		SxRNAIndividual.sequence_metadata_constraints = sequence_metadata_constrainst;
		SxRNAIndividual.motif_structure = motif_structure;
		SxRNAIndividual.begin_absolute_motif_struct_pos = SxRNAIndividual.motif_structure.indexOf('(');
		SxRNAIndividual.end_absolute_motif_struct_pos = SxRNAIndividual.motif_structure.lastIndexOf(')');
		
		//this is a hack to use current motif definition but allow for more complex multiple defs later
		SxRNAIndividual.motifs = new ArrayList<String>();
		SxRNAIndividual.motifs.add(SxRNAIndividual.motif_structure);
		//end hack
		
		int[] crossover_points = new int[2];
		crossover_points[0] = SxRNAIndividual.begin_absolute_motif_struct_pos;
		crossover_points[1] = SxRNAIndividual.end_absolute_motif_struct_pos;
				
		Individual.setDefaultCrossover_points(crossover_points);
	}



	@Override
	public Object clone() throws CloneNotSupportedException{
		SxRNAIndividual newIndividual = (SxRNAIndividual)super.clone();
		return newIndividual;
		
	}

	@Override
	public Individual getChildInstance() {
		return new SxRNAIndividual();
	}
	
	/**
	 * Returns a valid FASTA nucleotide sequence (sequence only) representing the current 
	 * set of alleles in the sxRNAIndividual's chromosome.
	 * @return the FASTA sequence the object represents.
	 */
	public String getFASTA_representation() {
		StringBuffer sb = new StringBuffer();
		for(int i =0; i<this.getGenes().size();i++){
			RibonucleotideGene currentGene = (RibonucleotideGene)this.getGenes().get(i); 
			if(!currentGene.toString().equals(" ")){
				sb.append(currentGene.toString());
				if(i==SxRNAIndividual.begin_absolute_motif_struct_pos){
					this.begin_motif_pos_in_fasta = sb.length()-1;					
				}else if(i==SxRNAIndividual.end_absolute_motif_struct_pos){
					this.end_motif_pos_in_fasta = sb.length()-1;
				}				
			}
		}
		return sb.toString();
	}
	/**
	 * Returns a full FASTA representation (with header line) representing the 
	 * sxRNAIndividual's ID and the set of alleles in its chromosome.
	 * @return the multi-line full FASTA representation the object represents   
	 */
	public String getFullFastaRepresentation(){
		StringBuffer sb = new StringBuffer();
		sb.append(">sxrna_");
		sb.append(this.getId());
		sb.append("\n");
		sb.append(this.getFASTA_representation());
		sb.append("\n");
		return sb.toString();
	}
	


	/**
	 * Converts the sequence_nucleotide_constraints FASTA code sequence into 
	 * permissible nucleotide arrays to create RibonucleotideGenes and determines
	 * from metadata what genes may have blank alleles.
	 */
	private void createNucleotideGenes() {
		ArrayList<Gene> newGenes = new ArrayList<Gene>();
		
			for(int i=0;i<SxRNAIndividual.sequence_nucleotide_constraints.length();i++){
				char currenFASTACode = SxRNAIndividual.sequence_nucleotide_constraints.charAt(i);
				char currentMetadata = SxRNAIndividual.sequence_metadata_constraints.charAt(i);
				RibonucleotideGene newGene = null;
				char[] permissableNucleotides = null;
				boolean mayBeBlank = false;
				
				switch(currenFASTACode){
					case 'A': permissableNucleotides = Ribonucleotides.a;
						break;
					case 'a': permissableNucleotides = Ribonucleotides.a;
						break;
					case 'C': permissableNucleotides = Ribonucleotides.c;
						break;
					case 'c': permissableNucleotides = Ribonucleotides.c;
						break;	
					case 'G': permissableNucleotides = Ribonucleotides.g;
						break;
					case 'g': permissableNucleotides = Ribonucleotides.g;
						break;
					case 'T': permissableNucleotides = Ribonucleotides.t;
						break;
					case 't': permissableNucleotides = Ribonucleotides.t;
						break;
					case 'U': permissableNucleotides = Ribonucleotides.u;
						break;
					case 'u': permissableNucleotides = Ribonucleotides.u;
						break;
					case 'R': permissableNucleotides = Ribonucleotides.r;
						break;
					case 'r': permissableNucleotides = Ribonucleotides.r;
						break;
					case 'Y': permissableNucleotides = Ribonucleotides.y;
						break;
					case 'y': permissableNucleotides = Ribonucleotides.y;
						break;
					case 'K': permissableNucleotides = Ribonucleotides.k;
						break;
					case 'k': permissableNucleotides = Ribonucleotides.k;
						break;
					case 'M': permissableNucleotides = Ribonucleotides.m;
						break;
					case 'm': permissableNucleotides = Ribonucleotides.m;
						break;
					case 'S': permissableNucleotides = Ribonucleotides.s;
						break;
					case 's': permissableNucleotides = Ribonucleotides.s;
						break;
					case 'W': permissableNucleotides = Ribonucleotides.w;
						break;
					case 'w': permissableNucleotides = Ribonucleotides.w;
						break;
					case 'B': permissableNucleotides = Ribonucleotides.b;
						break;
					case 'b': permissableNucleotides = Ribonucleotides.b;
						break;
					case 'D': permissableNucleotides = Ribonucleotides.d;
						break;
					case 'd': permissableNucleotides = Ribonucleotides.d;
						break;
					case 'H': permissableNucleotides = Ribonucleotides.h;
						break;
					case 'h': permissableNucleotides = Ribonucleotides.h;
						break;
					case 'V': permissableNucleotides = Ribonucleotides.v;
						break;
					case 'v': permissableNucleotides = Ribonucleotides.v;
						break;
					case 'N': permissableNucleotides = Ribonucleotides.n;
						break;
					case 'n': permissableNucleotides = Ribonucleotides.n;
						break;
					//this should not happen, perhaps an exception would be more appropriate
					default:  permissableNucleotides = Ribonucleotides.n;
						break;
				}
				switch(currentMetadata){
					case 'Y': mayBeBlank = true;
						break;
					case 'y': mayBeBlank = true;
						break;
					default: mayBeBlank = false;
						break;				
				}
				//System.out.println("mayBeBlank:"+mayBeBlank);
				newGene = new RibonucleotideGene(permissableNucleotides, mayBeBlank);
				newGenes.add(newGene);
			}
		
		
		this.setGenes(newGenes);
	}
	/**
	 * Returns the RNAFold dot parantheses MFE result for the candidate.
	 * @return the MFE fold as a dot parantheses string
	 */
	public String getRnaFold() {
		return rnaFold;
	}

	/**
	 * Sets the RNAFold dot parantheses MFE result
	 * @param rnaFold 
	 *        the MFE fold as a dot parantheses string
	 */
	public void setRnaFold(String rnaFold) {
		this.rnaFold = rnaFold;
	}

	/**
	 * Returns the RNAcofold dot parantheses MFE result for the candidate 
	 * with the primary trigger.
	 * @return the MFE cofold as a dot parantheses string
	 */
	public String getRnaCofold() {
		return rnaCofold;
	}

	/**
	 * Sets the RNAcofold dot parantheses MFE result for the candidate
	 * @param rnaCofold
	 *        the MFE cofold as a dot parantheses string
	 */
	public void setRnaCofold(String rnaCofold) {
		this.rnaCofold = rnaCofold;
	}


	/** 
	 * Provides cleanup services that might be needed after cloning and mutating 
	 * an sxRNAIndividual.
	 */
	@Override
	protected void postMutationHousekeeping() {		
		this.setRnaFold(null);
		this.setRnaCofold(null);
		
	}

	/**
	 * Returns the beginning genome position for the motif structure in the sxRNAIndividual 
	 * as defined by the motif_constraints argument
	 * @return the fixed beginning position of the desired motif structure in this 
	 * sxRNAIndividual
	 */
	public static int getBegin_absolute_motif_struct_pos() {
		return begin_absolute_motif_struct_pos;
	}

	/**
	 * Returns the end genome position for the motif structure in the sxRNAIndividual 
	 * as defined by the motif_constraints argument
	 * @return the fixed end position of the desired motif structure in this 
	 * sxRNAIndividual
	 */
	public static int getEnd_absolute_motif_struct_pos() {
		return end_absolute_motif_struct_pos;
	}

	/**
	 * Returns the beginning string position for the motif structure in the sxRNAIndividual 
	 * as it is located in FASTA sequence representation. This may differ from the genomic 
	 * position due to incorporation of blank alleles. 
	 * @return the beginning position of the desired motif structure in this 
	 * sxRNAIndividual's FASTA sequence
	 */
	public int getBegin_motif_pos_in_fasta() {
		return begin_motif_pos_in_fasta;
	}

	/**
	 * Returns the end string position for the motif structure in the sxRNAIndividual 
	 * as it is located in FASTA sequence representation. This may differ from the genomic 
	 * position due to incorporation of blank alleles. 
	 * @return the end position of the desired motif structure in this 
	 * sxRNAIndividual's FASTA sequence
	 */
	public int getEnd_motif_pos_in_fasta() {
		return end_motif_pos_in_fasta;
	}

	/**
	 * A variable length list of motif structures that may satisfy the desired functional motif.
	 * This is meant to allow for extensibility, currently only the first position is used.
	 * @return the list of motif structures
	 */
	public static List<String> getMotifs() {
		return motifs;
	}

	/**
	 * Checks if the trigger RNA is predicted to hybridize acroos base of the motif region
	 * in the sxRNAIndividual.
	 * @return {@code true} if the trigger hybridizes across base of motif in the RNAcofold 
	 * MFE prediction, {@code false} otherwise
	 */
	public boolean trigger_spans_motif() {
		return trigger_spans_motif;
	}

	/**
	 * Sets variable telling whether this candidates cofold MFE shows trigger hybridized to both sides of 
	 * the motif region
	 * @param trigger_spans_motif
	 * 		the value telling whether this candidate's cofold MFE shows trigger binding across motif region
	 */
	public void setTrigger_spans_motif(boolean trigger_spans_motif) {
		this.trigger_spans_motif = trigger_spans_motif;
	}

	/**
	 * Returns the non target ncRNA's that are predicted to allow motif to form 
	 * @return the String list of non target ncRNAs that facilitate motif forming  
	 */
	public String getNonTargetPosHits() {
		return nonTargetPosHits;
	}

	/**
	 * Sets the String list of non targeted ncRNAs that seem to facilitate motif formation
	 * @param nonTargetPosHits
	 * 	the String list of non target ncRNAs that facilitate motif forming
	 */
	public void setNonTargetPosHits(String nonTargetPosHits) {
		this.nonTargetPosHits = nonTargetPosHits;
	}

	
	public String getNonTargetSuppressorHits() {
		return nonTargetSuppressorHits;
	}


	public void setNonTargetSuppressorHits(String nonTargetSuppressorHits) {
		this.nonTargetSuppressorHits = nonTargetSuppressorHits;
	}


	public String getNegTargetPosHits() {
		return negTargetPosHits;
	}


	public void setNegTargetPosHits(String negTargetPosHits) {
		this.negTargetPosHits = negTargetPosHits;
	}


	public String getNegTargetSuppressorHits() {
		return negTargetSuppressorHits;
	}


	public void setNegTargetSuppressorHits(String negTargetSuppressorHits) {
		this.negTargetSuppressorHits = negTargetSuppressorHits;
	}


	public String getPosTargetPosHits() {
		return posTargetPosHits;
	}


	public void setPosTargetPosHits(String posTargetPosHits) {
		this.posTargetPosHits = posTargetPosHits;
	}


	public String getPosTargetSuppressorHits() {
		return posTargetSuppressorHits;
	}


	public void setPosTargetSuppressorHits(String posTargetSuppressorHits) {
		this.posTargetSuppressorHits = posTargetSuppressorHits;
	}

	public double getPosTargetPosHitScore() {
		return posTargetPosHitScore;
	}


	public void setPosTargetPosHitScore(double posTargetPosHitScore) {
		this.posTargetPosHitScore = posTargetPosHitScore;
	}


	public double getPosTargetSuppressorHitScore() {
		return posTargetSuppressorHitScore;
	}


	public void setPosTargetSuppressorHitScore(double posTargetSuppressorHitScore) {
		this.posTargetSuppressorHitScore = posTargetSuppressorHitScore;
	}


	public double getNonTargetPosHitScore() {
		return nonTargetPosHitScore;
	}


	public void setNonTargetPosHitScore(double nonTargetPosHitScore) {
		this.nonTargetPosHitScore = nonTargetPosHitScore;
	}


	public double getNonTargetSuppressorHitScore() {
		return nonTargetSuppressorHitScore;
	}


	public void setNonTargetSuppressorHitScore(double nonTargetSuppressorHitScore) {
		this.nonTargetSuppressorHitScore = nonTargetSuppressorHitScore;
	}


	public double getNegTargetPosHitScore() {
		return negTargetPosHitScore;
	}


	public void setNegTargetPosHitScore(double negTargetPosHitScore) {
		this.negTargetPosHitScore = negTargetPosHitScore;
	}


	public double getNegTargetSuppressorHitScore() {
		return negTargetSuppressorHitScore;
	}


	public void setNegTargetSuppressorHitScore(double negTargetSuppressorHitScore) {
		this.negTargetSuppressorHitScore = negTargetSuppressorHitScore;
	}


	@Override
	public void setGenesByStringForTesting(String string_representation) {
		if(string_representation != null){
			for(int i=0;i<string_representation.length();i++){
				RibonucleotideGene current_gene = (RibonucleotideGene)this.getGenes().get(i);
				char current_char = string_representation.charAt(i);
				if(current_char==' '){
					current_gene.setCurrentAlleleForTesting(' ');	
				}else if(current_char=='a'||current_char=='A'){
					current_gene.setCurrentAlleleForTesting(Ribonucleotides.a[0]);
				}else if(current_char=='c'||current_char=='C'){
					current_gene.setCurrentAlleleForTesting(Ribonucleotides.c[0]);
				}else if(current_char=='g'||current_char=='G'){
					current_gene.setCurrentAlleleForTesting(Ribonucleotides.g[0]);
				}else if(current_char=='u'||current_char=='U'){
					current_gene.setCurrentAlleleForTesting(Ribonucleotides.u[0]);
				}else if(current_char=='t'||current_char=='T'){
					current_gene.setCurrentAlleleForTesting(Ribonucleotides.u[0]);
				}else{
					System.err.println("Unrecognized ribonucleotide string code for testing!["+current_char+"]");
					
					System.exit(1);
				}				
			}
		}
		
	}
	
	
	@Override
	public String toString() {
		StringBuffer sb = new StringBuffer();
		sb.append(super.toString());
		sb.append("\n ");
		sb.append("[SxRNAIndividual] - fold: "+this.rnaFold);
		sb.append("   cofold: "+ this.rnaCofold);
		sb.append("   nonTargetPosScore["+this.nonTargetPosHitScore+"]:" + this.nonTargetPosHits);
		sb.append("   nonTargetSuppressorScore["+this.nonTargetSuppressorHitScore+"]:" + this.nonTargetSuppressorHits);
		sb.append("   negTargetPosScore["+this.negTargetPosHitScore+"]:" + this.negTargetPosHits);
		sb.append("   negTargetSuppressorScore["+this.negTargetSuppressorHitScore+"]:" + this.negTargetSuppressorHits);
		sb.append("   posTargetPosScore["+this.posTargetPosHitScore+"]:" + this.posTargetPosHits);
		sb.append("   posTargetSuppressorScore["+this.posTargetSuppressorHitScore+"]:" + this.posTargetSuppressorHits);
		return sb.toString();
	}



	/**
	 * @return
	 */
	public boolean isMotif_present_in_solo_fold() {
		return motif_present_in_solo_fold;
	}


	/**
	 * @param motif_present_in_solo_fold
	 */
	public void setMotif_present_in_solo_fold(boolean motif_present_in_solo_fold) {
		this.motif_present_in_solo_fold = motif_present_in_solo_fold;
	}


	/**
	 * @return
	 */
	public boolean isMotif_present_in_cofold() {
		return motif_present_in_cofold;
	}


	/**
	 * @param motif_present_in_cofold
	 */
	public void setMotif_present_in_cofold(boolean motif_present_in_cofold) {
		this.motif_present_in_cofold = motif_present_in_cofold;
	}
	

	/**
	 * 
	 * @return
	 */
	public String getSupplementaryCensusData() {
		StringBuffer sb = new StringBuffer();
		sb.append("|nonTargetPosScore["+this.nonTargetPosHitScore+"]");
		sb.append("|nonTargetSuppressScore["+this.nonTargetSuppressorHitScore+"]");
		sb.append("|negTargetPosScore["+this.negTargetPosHitScore+"]");
		sb.append("|negTargetSuppressScore["+this.negTargetSuppressorHitScore+"]");
		sb.append("|posTargetPosScore["+this.posTargetPosHitScore+"]");
		sb.append("|posTargetSuppressScore["+this.posTargetSuppressorHitScore+"]");
		return sb.toString();
	}







	
	
	
}
