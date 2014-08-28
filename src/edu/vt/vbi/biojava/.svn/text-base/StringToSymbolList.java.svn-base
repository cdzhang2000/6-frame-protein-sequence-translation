package edu.vt.vbi.biojava;

import org.biojava.bio.seq.*;
import org.biojava.bio.symbol.*;

public class StringToSymbolList {
	
	public void stringToSymbolList(){
		

		try {
			// create a DNA SymbolList from a String
			SymbolList dna = DNATools.createDNA("atcggtcggctta");
			
			int s=dna.length();
			System.out.println("DNA sublist="+ dna.subStr(1, s));
			
			int c=1;
			
			while(s>0 && c++<2){
				System.out.print("Alphabet ="+ dna.getAlphabet().getName());
				System.out.print("\tSymbol ="+ dna.symbolAt(s).getName());
				System.out.println();
			}

			c=1;
			// create a RNA SymbolList from a String
			SymbolList rna = RNATools.createRNA("auugccuacauaggc");
			
			s=rna.length();
			System.out.println("RNA sublist="+ rna.subStr(1, s));
			while(s>0&& c++<2){
				System.out.print("Alphabet ="+ rna.getAlphabet().getName());
				System.out.print("\tSymbol ="+ rna.symbolAt(s).getName());
				System.out.println();
			}
			c=1;
			// create a Protein SymbolList from a String
			SymbolList aa = ProteinTools.createProtein("AGFAVENDSA");
			s=aa.length();
			while(s>0&& c++<2){
				System.out.print("Alphabet ="+ aa.getAlphabet().getName());
				System.out.print("\tSymbol ="+ aa.symbolAt(s).getName());
				System.out.println();
			}

			
		} catch (IllegalSymbolException ex) {
			ex.printStackTrace();
		}
	}
	
	public static void main(String[] args) {
		
		StringToSymbolList s= new StringToSymbolList();
		s.stringToSymbolList();
	}
}
