package edu.vt.vbi.biojava;
import java.io.*;
import java.util.*;
 
import org.biojava.bio.*;
import org.biojava.bio.seq.*;
import org.biojava.bio.seq.io.*;

public class ListAnnotations {
	public static void main(String[] args) {
		 
	    try {
	      //read in an EMBL Record
	      BufferedReader br = new  BufferedReader(new FileReader(args[0]));
	 
	      //for each sequence list the annotations
	      for(SequenceIterator seqs = SeqIOTools.readEmbl(br); seqs.hasNext(); ){
	        Annotation anno = seqs.nextSequence().getAnnotation();
	 
	        //print each key value pair
	        for (Iterator i = anno.keys().iterator(); i.hasNext(); ) {
	          Object key = i.next();
	          System.out.println(key +" : "+ anno.getProperty(key));
	        }
	      }
	    }
	    catch (Exception ex) {
	      ex.printStackTrace();
	    }
	  }

}
