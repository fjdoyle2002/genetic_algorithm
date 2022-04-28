package io.github.fjdoyle2002.util;

/**
 * @author fd299212
 * Class represents a non coding RNA 
 */
public class NonCodingRNA extends Object{
	private String id = null;
	private String header = null;
	private String sequence = null;	
	
	public NonCodingRNA(String header, String sequence) {
		super();
		this.setHeader(header);
		if(!(null==header)){
			int cutIndex = header.indexOf(" ");
			if(!(cutIndex>0)){
				cutIndex = header.length();
			}
			this.setId(header.substring(0,cutIndex));
		}
		this.setSequence(sequence);
	}
	/**
	 * @return the id
	 */
	public String getId() {
		return id;
	}
	private void setId(String id){
		this.id= id;
	}
	/**
	 * @return the header
	 */
	public String getHeader() {
		return header;
	}
	/**
	 * @param header the header to set
	 */
	private void setHeader(String header) {
		this.header = header;
	}
	/**
	 * @return the sequence
	 */
	public String getSequence() {
		return sequence;
	}
	/**
	 * @param sequence the sequence to set
	 */
	private void setSequence(String sequence) {
		this.sequence = sequence;
	}
	/**
	 * 
	 * @return
	 */
	public String getFullFastaRepresentation(){
		StringBuffer sb = new StringBuffer();
		sb.append(">");
		sb.append(this.header);
		sb.append("\n");
		sb.append(this.sequence);
		sb.append("\n");
		return sb.toString();
	}
	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		return "NonCodingRNA [header=" + header + ", sequence=" + sequence
				+ "]";
	}

	
	
}
