package edu.vt.vbi.biojava;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.NoSuchElementException;

import org.biojava.bio.Annotation;
import org.biojava.bio.BioException;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.RNATools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.SequenceTools;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.SymbolList;
import org.biojavax.bio.seq.RichSequence;
/**
 * 
 * @author CZHANG
 * Take FASTA file and translate DNA sequence into 6-frame proteins sequence. 
 * Protein sequence is written into a file which can be named in the main method.    
 *
 */

public class HexTranslate2 {

	// private Vector<ORF> orfs = null;

	public HexTranslate2() {
		// orfs = new Vector<ORF>();
	}

	public void myTranslate(String fileName, String type, String output) {
		try {
			SymbolTokenization toke = AlphabetManager.alphabetForName(type)
					.getTokenization("token");
			BufferedReader br = new BufferedReader(new FileReader(fileName));
			SequenceIterator seqi = RichSequence.IOTools.readFasta(br, toke,
					null);

			String seqName="";

			File file = new File(output);

			BufferedWriter writer = new BufferedWriter(new FileWriter(file));

			// for each sequence
			if (seqi.hasNext()) { //change if to while, it can read multiple seq in one file

				Sequence seq = seqi.nextSequence();
				
				seqName = seq.getName();
				
				System.out.println("sequence name ="+seqName);
				
				// for each frame
				for (int i = 0; i < 3; i++) {

					forwardTranslate(seq, i, seqName, writer);

					backwardTranslate(seq, i, seqName, writer);

				}
			}
			br.close();
			writer.flush();
			writer.close();
			
			System.out.println("END name="+seqName);
			
			
		} catch (IOException e) {
			e.printStackTrace();
		} catch (IllegalAlphabetException e) {
			e.printStackTrace();
		} catch (NoSuchElementException e) {
			e.printStackTrace();
		} catch (BioException e) {
			e.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	public void forwardTranslate(Sequence seq, int i, String seqName,
			BufferedWriter w) throws IllegalAlphabetException,
			NoSuchElementException, BioException, IOException {
		SymbolList prot;
		Sequence trans;

		int length = seq.length() - (seq.length() - i) % 3;

		SymbolList syms = seq.subList(i + 1, length);

		// if it is DNA transcribe it to RNA
		if (syms.getAlphabet() == DNATools.getDNA()) {
			syms = DNATools.toRNA(syms);
		}

		// output forward translation to STDOUT
		prot = RNATools.translate(syms);

		trans = SequenceTools.createSequence(prot, "", seqName
				+ "TranslationFrame: +" + (i + 1), Annotation.EMPTY_ANNOTATION);

		subStr(trans, i + 1, "forward", seqName, w);

	}

	public void backwardTranslate(Sequence seq, int i, String seqName,
			BufferedWriter w) throws Exception {
		SymbolList prot;
		Sequence trans;

		int seqlength = seq.length();

		SymbolList syms = seq.subList(1, seqlength);

		// if it is DNA transcribe it to RNA
		if (syms.getAlphabet() == DNATools.getDNA()) {
			syms = DNATools.toRNA(syms);
		}

		syms = RNATools.reverseComplement(syms);

		int rlength = syms.length() - (syms.length() - i) % 3;

		syms = syms.subList(i + 1, rlength);

		prot = RNATools.translate(syms);

		trans = SequenceTools.createSequence(prot, "", seqName
				+ "TranslationFrame: +" + (i + 1), Annotation.EMPTY_ANNOTATION);

		int frame = i + 1;
		System.out.println("frame= " + frame + "\tDNA length=" + seqlength
				+ "\t rna length" + rlength);

		subStrBackward(trans, frame, seqName, rlength, seqlength, w);

	}

	public static void main(String[] args) {

		HexTranslate2 hex = new HexTranslate2();
		// String fileName = "H:/fasta/test.fna";
		//String fileName = "H:/fasta/NC_009089.fna";
		//String output = "H:/fasta/NC_009089_6frame.fna";
		String fileName = "H:/fasta/difficile_R20291.fna";
		String output = "H:/fasta/difficile_R20291_6frame.fna";
		String seqType = "DNA";
		hex.myTranslate(fileName, seqType, output);
	}

	public void subStr(Sequence seq, int frame, String strand, String seqName,
			BufferedWriter w) throws IOException {

		int length = seq.length();
		int start = frame;
		int stop = frame;

		int order = 1;

		System.out.println("seq length= " + length);

		String temp = seq.subStr(1, length);

		// set start position
		int s = 0;

		char t = '|';

		StringBuffer sb = new StringBuffer();

		while (s < length) {

			t = temp.charAt(s);

			// reset stop position
			stop = (s + 1) * 3 + frame - 1;

			if (t == '*') {

				ORF orf = new ORF();

				orf.setStart(start);

				// reset start
				start = stop + 1;

				// minus stop code position
				orf.setStop(stop - 3);

				orf.setFrame(strand + "_" + frame);

				orf.setId(seqName + "_" + orf.getStart() + "_" + orf.getStop()
						+ "_" + orf.getFrame());

				orf.setAaseq(sb.toString());

				orf.setOrder(order);
				orf.setS(s + 1);
				order = order + 1;
				orf.display();

				// orfs.add(orf);
				if (orf.getAaseq().length() > 30) {
					orf.writeToFile(w);
				}
				sb = null;
				sb = new StringBuffer();

			} else if (s == length - 1) {
				ORF orf = new ORF();
				orf.setStart(start);
				orf.setStop(stop);
				orf.setFrame(strand + "_" + frame);

				orf.setId(seqName + "_" + orf.getStart() + "_" + orf.getStop()
						+ "_" + orf.getFrame());
				sb.append(t);
				orf.setAaseq(sb.toString());
				orf.setOrder(order);
				orf.setS(s + 1);
				orf.display();

				// orfs.add(orf);
				if (orf.getAaseq().length() > 30) {
					orf.writeToFile(w);
				}
				// orfs.add(orf);
			} else {
				sb.append(t);
			}
			s = s + 1;
		}

	}

	public void subStrBackward(Sequence seq, int frame, String seqName,
			int rlength, int slength, BufferedWriter w) throws IOException {

		// protein sequence length
		int length = seq.length();

		// DNA sequence length
		int start = slength - (frame - 1);

		// set initial stop position
		int stop = rlength;

		int order = 1;

		// put protein sequence into a String
		String temp = seq.subStr(1, length);

		// String starts position 0
		int s = 0;

		char t = '|';

		StringBuffer sb = new StringBuffer();

		while (s < length) {

			t = temp.charAt(s);

			stop = s * 3;

			if (t == '*') {

				ORF orf = new ORF();
				orf.setStart(start);

				// end position
				orf.setStop(slength - (stop - 1 + (frame - 1)));

				// reset start position
				start = slength - (stop + 3 + (frame - 1));

				orf.setFrame("backward_" + frame);

				orf.setId(seqName + "_" + orf.getStart() + "_" + orf.getStop()
						+ "_" + orf.getFrame());

				// put sequence in the StringBuffer
				orf.setAaseq(sb.toString());

				orf.setOrder(order);
				orf.setS(s);
				order = order + 1;
				orf.display();
				// orfs.add(orf);

				// orfs.add(orf);
				if (orf.getAaseq().length() > 30) {
					orf.writeToFile(w);
				}
				sb = null;
				sb = new StringBuffer();

			} else if (s == length - 1) {

				sb.append(t);

				ORF orf = new ORF();

				orf.setStart(start);

				orf.setStop(slength - (stop + 3 + (frame - 1) + 1));

				orf.setFrame("backward_" + frame);

				orf.setId(seqName + "_" + orf.getStart() + "_" + orf.getStop()
						+ "_" + orf.getFrame());

				orf.setAaseq(sb.toString());

				orf.setOrder(order);

				orf.display();

				// orfs.add(orf);
				if (orf.getAaseq().length() > 30) {
					orf.writeToFile(w);
				}
				// orfs.add(orf);

			} else {
				sb.append(t);
			}
			// always move forward 1 position
			s = s + 1;
		}

	}

}
