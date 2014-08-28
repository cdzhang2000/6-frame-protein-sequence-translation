package edu.vt.vbi.biojava;

import java.io.BufferedReader;
import java.io.FileReader;
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

public class HexTranslate {

	//private Vector<ORF> orfs = null;

	public HexTranslate() {
	
	}

	public void myTranslate(String fileName, String type) {
		try {
			SymbolTokenization toke = AlphabetManager.alphabetForName(type)
					.getTokenization("token");
			BufferedReader br = new BufferedReader(new FileReader(fileName));
			SequenceIterator seqi = RichSequence.IOTools.readFasta(br, toke,
					null);

			String seqName;

			// for each sequence
			while (seqi.hasNext()) {

				Sequence seq = seqi.nextSequence();
				seqName = seq.getName();

				// for each frame
				for (int i = 0; i < 3; i++) {

					forwardTranslate(seq, i, seqName);

					backwardTranslate(seq, i, seqName);

				}
			}
			br.close();

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

	public void backwardTranslate(Sequence seq, int i, String seqName)
			throws Exception {
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

		subStrBackward(trans, frame, seqName, rlength, seqlength);

	}

	public void forwardTranslate(Sequence seq, int i, String seqName)
			throws IllegalAlphabetException, NoSuchElementException,
			BioException {
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

		subStr(trans, i + 1, "forward", seqName);

	}

	public static void main(String[] args) {

		HexTranslate hex = new HexTranslate();
		String fileName = "H:/fasta/NC_009089.fna";
		String seqType = "DNA";
		hex.myTranslate(fileName, seqType);
	}

	public void subStr(Sequence seq, int frame, String strand, String seqName) {

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
				//orfs.add(orf);
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
			} else {
				sb.append(t);
			}
			s = s + 1;
		}

	}

	public void subStrBackward(Sequence seq, int frame, String seqName,
			int rlength, int slength) {

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

				orf.setFrame("backword_" + frame);

				orf.setId(seqName + "_" + orf.getStart() + "_" + orf.getStop()
						+ "_" + orf.getFrame());

				// put sequence in the StringBuffer
				orf.setAaseq(sb.toString());

				orf.setOrder(order);
				orf.setS(s);
				order = order + 1;
				orf.display();
				//orfs.add(orf);
				sb = null;
				sb = new StringBuffer();

			} else if (s == length - 1) {

				sb.append(t);

				ORF orf = new ORF();

				orf.setStart(start);

				orf.setStop(slength - (stop + 3 + (frame - 1) + 1));

				orf.setFrame("backword_" + frame);

				orf.setId(seqName + "_" + orf.getStart() + "_" + orf.getStop()
						+ "_" + orf.getFrame());

				orf.setAaseq(sb.toString());

				orf.setOrder(order);

				orf.display();

				//orfs.add(orf);

			} else {
				sb.append(t);
			}
			// always move forward 1 position
			s = s + 1;
		}

	}

}
