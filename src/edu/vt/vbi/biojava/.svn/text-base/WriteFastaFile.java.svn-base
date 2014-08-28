package edu.vt.vbi.biojava;

import java.io.IOException;

import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.db.HashSequenceDB;
import org.biojava.bio.seq.db.SequenceDB;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojavax.bio.seq.RichSequence;

public class WriteFastaFile {
	
	private static void printSequenceDB() {
		
		SequenceDB db = new HashSequenceDB();
		
		Sequence dna1;
		Sequence dna2;
		try {
			dna1 = DNATools.createDNASequence("atgctgtgg", "dna_1");
			dna2 = DNATools.createDNASequence("atgctgctt", "dna_2");
			db.addSequence(dna1);
		        db.addSequence(dna2);
			RichSequence.IOTools.writeFasta(System.out, db.sequenceIterator(), null);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	 /*
     * SeqIOTools also has a method that takes a single sequence so you don't
     * have to make a SequenceDB
     */
    	private static void printSingleSequence(){
		Sequence dna;
		try {
			dna = DNATools.createDNASequence("atgctg", "dna_1");
			RichSequence.IOTools.writeFasta(System.out, dna, null);
		} catch (IllegalSymbolException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}		
	}
}
