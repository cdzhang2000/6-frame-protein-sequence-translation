package edu.vt.vbi.biojava;
/**
 * Class to load an EMBL or Genbank sequence file and output annotation information about the CDS features.
 */
 
//Java libraries
import java.io.*;
import java.util.*;
//BioJava libraries
import org.biojava.bio.*;
import org.biojava.bio.seq.*;
import org.biojava.bio.seq.io.*;
//BioJava extension libraries
import org.biojavax.*;
import org.biojavax.ontology.*;
import org.biojavax.bio.*;
import org.biojavax.bio.seq.*;


public class ExtractInformation {
	 //Create the RichSequence object
	  RichSequence richSeq;
	 
	  //ExtractInformation constructor
	  public ExtractInformation(String fileName){
	    //Load the sequence file
	    try {
	      richSeq = RichSequence.IOTools.readGenbankDNA(new BufferedReader(new FileReader(fileName)),null).nextRichSequence();
	    }
	    catch(FileNotFoundException fnfe){
	      System.out.println("FileNotFoundException: " + fnfe);
	    }
	    catch(BioException bioe1){
	      System.err.println("Not a Genbank sequence trying EMBL");
	      try  {
	        richSeq = RichSequence.IOTools.readEMBLDNA(new BufferedReader(new FileReader(fileName)),null).nextRichSequence();
	      }
	      catch(BioException bioe2){
	        System.err.println("Not an EMBL sequence either");
	        System.exit(1);
	      }
	      catch(FileNotFoundException fnfe){
	        System.out.println("FileNotFoundException: " + fnfe);
	      }
	    }
	    //Filter the sequence on CDS features
	    FeatureFilter ff = new FeatureFilter.ByType("CDS");
	    FeatureHolder fh = richSeq.filter(ff);
	 
	    //Iterate through the CDS features
	    for (Iterator<Feature> i = fh.features(); i.hasNext();){
	      RichFeature rf = (RichFeature) i.next();
	 
	      //Get the strand orientation of the feature
	      char featureStrand = rf.getStrand().getToken();
	 
	      //Get the location of the feature
	      String featureLocation = rf.getLocation().toString();
	 
	      //Get the annotation of the feature
	      RichAnnotation ra = (RichAnnotation)rf.getAnnotation();
	 
	      //Use BioJava defined ComparableTerms 
	      ComparableTerm geneTerm = new RichSequence.Terms().getGeneNameTerm();
	      ComparableTerm synonymTerm = new RichSequence.Terms().getGeneSynonymTerm();
	      //Create the required additional ComparableTerms
	      ComparableTerm locusTerm = RichObjectFactory.getDefaultOntology().getOrCreateTerm("locus_tag");
	      ComparableTerm productTerm = RichObjectFactory.getDefaultOntology().getOrCreateTerm("product");
	      ComparableTerm proteinIDTerm = RichObjectFactory.getDefaultOntology().getOrCreateTerm("protein_id");
	 
	      //Create empty strings
	      String gene = "";
	      String locus = "";
	      String product = "";
	      String geneSynonym = "";
	      String proteinID = "";
	 
	      //Iterate through the notes in the annotation 
	      for (Iterator <Note> it = ra.getNoteSet().iterator(); it.hasNext();){
	        Note note = it.next();
	 
	      //Check each note to see if it matches one of the required ComparableTerms
	        if(note.getTerm().equals(locusTerm)){
	          locus = note.getValue().toString();
	        }
	        if(note.getTerm().equals(productTerm)){
	          product = note.getValue().toString();
	        }
	        if(note.getTerm().equals(geneTerm)){
	          gene = note.getValue().toString();
	        }
	        if(note.getTerm().equals(synonymTerm)){
	          geneSynonym = note.getValue().toString();
	        }
	        if(note.getTerm().equals(proteinIDTerm)){
	          proteinID = note.getValue().toString();
	        }
	      }
	      //Outout the feature information
	      System.out.println(locus + "  " + gene + "  " + geneSynonym + "  " + proteinID + "  " + product + "  " + featureStrand + "  " + featureLocation);
	    }
	  }
	 
	  //Main method
	  public static void main(String args []){
	    if (args.length != 1){
	      System.out.println("Usage: java ExtractInformation <file in Genbank or EMBL format>");
	      System.exit(1);
	    }
	    else {
	      new ExtractInformation(args[0]);
	    }
	  }
}
