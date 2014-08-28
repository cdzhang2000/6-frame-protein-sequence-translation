package edu.vt.vbi.biojava;
import org.biojava.bio.seq.io.*;
import org.biojava.bio.seq.*;
import org.biojavax.Namespace;
import org.biojavax.RichObjectFactory;
import org.biojavax.SimpleNamespace;
import org.biojavax.bio.seq.RichSequence;
import org.biojavax.bio.seq.RichSequenceIterator;
import org.biojavax.bio.seq.io.RichSequenceBuilderFactory;
import org.biojavax.bio.seq.io.RichSequenceFormat;
import org.biojavax.bio.seq.io.RichStreamReader;
import org.biojavax.bio.seq.io.RichStreamWriter;

import java.io.*;

public class WriteToFasta {
	/**
	   * This program will read any file supported by SeqIOTools it takes two
	   * arguments, the first is the file name the second is the int constant
	   * for the file type in SeqIOTools. See SeqIOTools for possible file types.
	   * The constants used are:
	   * UNKNOWN = 0;
	   * FASTADNA = 1;
	   * FASTAPROTEIN = 2;
	   * EMBL = 3;
	   * GENBANK = 4;
	   * SWISSPROT = 5;
	   * GENPEPT = 6;
	   * MSFDNA = 7;
	   * FASTAALIGNDNA = 9;
	   * MSFPROTEIN = 10;
	   * FASTAALIGNPROTEIN = 11;
	   * MSF = 12;
	   *
	   */
	  public static void main(String[] args) {
	    try {
	      //prepare a BufferedReader for file io
	      BufferedReader br = new BufferedReader(new FileReader(args[0]));
	 
	      //get the int constant for the file type
	      int fileType = Integer.parseInt(args[1]);
	 
	      /*
	       * get a Sequence Iterator over all the sequences in the file.
	       * SeqIOTools.fileToBiojava() returns an Object. If the file read
	       * is an alignment format like MSF and Alignment object is returned
	       * otherwise a SequenceIterator is returned.
	       */
	      SequenceIterator iter =
	          (SequenceIterator)SeqIOTools.fileToBiojava(fileType, br);
	 
	      //and now write it all to FASTA, (you can write to any OutputStream, not just System.out)
	      SeqIOTools.writeFasta(System.out, iter);
	    }
	    catch (Exception ex) {
	      ex.printStackTrace();
	    }
	  }
	  
	  public void writeToEMBL()throws Exception{
		// an input GenBank file
		  BufferedReader br = new BufferedReader(new FileReader("myGenbank.gbk"));  
		  // a namespace to override that in the file
		  Namespace ns = RichObjectFactory.getDefaultNamespace();                   
		  // we are reading DNA sequences
		  RichSequenceIterator seqs = RichSequence.IOTools.readGenbankDNA(br,ns);   
		  // write the whole lot in EMBLxml format to standard out
		  RichSequence.IOTools.writeEMBLxml(System.out, seqs, ns);
	  }
	  
	  /**
	   * read some DNA sequences from a GenBank file and write them out to standard output (screen) as FASTA
	   * @throws Exception
	   */
	  public void writeToFasta()throws Exception {
		// sequences will be DNA sequences
		  SymbolTokenization dna = DNATools.getDNA().getTokenization("token");        
		  // read Genbank
		  RichSequenceFormat genbank = (RichSequenceFormat) new GenbankFormat();                           
		  // write FASTA
		  RichSequenceFormat fasta = (RichSequenceFormat) new FastaFormat();                               
		  // compress only longer sequences
		  RichSequenceBuilderFactory factory = RichSequenceBuilderFactory.THRESHOLD;  
		  // read/write everything using the 'bloggs' namespace
		  Namespace bloggsNS = (Namespace) RichObjectFactory.getObject(
		                          SimpleNamespace.class, 
		                          new Object[]{"bloggs"} 
		                       );                                                     
		   
		  // read seqs from "mygenbank.file"
		  BufferedReader input = new BufferedReader(new FileReader("mygenbank.file"));
		  // write seqs to STDOUT
		  OutputStream output = System.out;                                           
		   
		  RichStreamReader seqsIn = new RichStreamReader(input,genbank,dna,factory,bloggsNS);
		  RichStreamWriter seqsOut = new RichStreamWriter(output,fasta);
		  // one-step Genbank to Fasta conversion!
		  seqsOut.writeStream(seqsIn,bloggsNS);
	  }
	  
	  

}
