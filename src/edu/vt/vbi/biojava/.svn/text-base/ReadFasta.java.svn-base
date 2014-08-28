package edu.vt.vbi.biojava;

import java.io.*;
import java.util.*;
 
import org.biojava.bio.*;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.db.*;
import org.biojava.bio.seq.io.*;
import org.biojava.bio.symbol.*;
import org.biojavax.Namespace;
import org.biojavax.RichObjectFactory;
import org.biojavax.bio.seq.RichSequence;
import org.biojavax.bio.seq.RichSequenceIterator;


public class ReadFasta {
	 /**
	   * The program takes two args: the first is the file name of the Fasta file.
	   * The second is the name of the Alphabet. Acceptable names are DNA RNA or PROTEIN.
	 * @throws IOException 
	 * @throws BioException 
	   */
	  public void readFasta(String fileName, String seqType) throws IOException, BioException {
	 
	    try {
	    
	   // an input GenBank file
	      BufferedReader br = new BufferedReader(new FileReader(fileName));  
	      // a namespace to override that in the file
	      Namespace ns = RichObjectFactory.getDefaultNamespace();                   
	      // we are reading DNA sequences
	      RichSequenceIterator seqs = RichSequence.IOTools.readGenbankDNA(br,ns);   
	      // write the whole lot in EMBLxml format to standard out
	      RichSequence.IOTools.writeEMBLxml(System.out, seqs, ns);
	     
	      while (seqs.hasNext()) {
	          RichSequence rs = seqs.nextRichSequence();
	          // write it in EMBL format to standard out
	          RichSequence.IOTools.writeEMBL(System.out, rs, ns);                   
	      }
	      
	    }
	    catch (NoSuchElementException ex) {
	      //no fasta sequences in the file
	      ex.printStackTrace();
	    }catch (FileNotFoundException ex) {
	      //problem reading file
	      ex.printStackTrace();
	    }
	  }
	  
	  /**
	   * This program will read any file supported by SeqIOTools it takes three
	   * arguments, the first is the file name the second is the name of
	   * a file format supported by SeqIOTools. eg fasta, genbank etc.
	   * The third argument is the alphabet (eg dna, rna, protein).
	   *
	   * Both the format and alphabet names are case insensitive.
	   *
	   */
	  public void readFasta(String fileName, String source, String seqType) {
	    try {
	      //prepare a BufferedReader for file io
	      BufferedReader br = new BufferedReader(new FileReader(fileName));
	 
	      String format = source;
	      String alphabet = seqType;
	 
	      /*
	       * get a Sequence Iterator over all the sequences in the file.
	       * SeqIOTools.fileToBiojava() returns an Object. If the file read
	       * is an alignment format like MSF and Alignment object is returned
	       * otherwise a SequenceIterator is returned.
	       */
	      SequenceIterator iter =
	          (SequenceIterator)SeqIOTools.fileToBiojava(format, alphabet, br);
	    }
	    catch (FileNotFoundException ex) {
	      //can't find file specified by args[0]
	      ex.printStackTrace();
	    }catch (BioException ex) {
	      //error parsing requested format
	      ex.printStackTrace();
	    }
	  }

	  
}
