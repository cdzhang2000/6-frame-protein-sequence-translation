package edu.vt.vbi.biojava;

import org.biojava.bio.seq.*;
import org.biojava.bio.symbol.*;

public class StringToSequence {

	public void toseq() {
		try {
			// create a DNA sequence with the name dna_1
			Sequence dna = DNATools.createDNASequence("atgctg", "dna_1");

			// create an RNA sequence with the name rna_1
			Sequence rna = RNATools.createRNASequence("augcug", "rna_1");

			// create a Protein sequence with the name prot_1
			Sequence prot = ProteinTools
					.createProteinSequence("AFHS", "prot_1");
		} catch (IllegalSymbolException ex) {
			// an exception is thrown if you use a non IUB symbol
			ex.printStackTrace();
		}
	}
	
	public static void main(String args[]){
		
		StringToSequence sts=new StringToSequence();
		sts.toseq();
		
	}
}
