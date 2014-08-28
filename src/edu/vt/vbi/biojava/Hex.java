package edu.vt.vbi.biojava;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.io.ObjectInputStream.GetField;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.Vector;

import org.biojava.bio.Annotation;
import org.biojava.bio.BioException;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.RNATools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.SequenceTools;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.biojavax.bio.seq.RichSequence;

public class Hex {

	private Vector<ORF> orfs = null;

	public Hex() {
		orfs = new Vector<ORF>();
	}

	public void readFasta(String fileName, String seqType) {
		String filename = fileName;
		String type = seqType;
		try {
			SymbolTokenization toke = AlphabetManager.alphabetForName(type)
					.getTokenization("token");

			BufferedReader br = new BufferedReader(new FileReader(filename));

			// SequenceIterator seqi = RichSequence.IOTools.readFasta(br, toke,
			// null);
			SequenceIterator seqi = RichSequence.IOTools.readFastaDNA(br, null);

			String vbi_id = "";

			int c = 1;

			// for each sequence
			while (seqi.hasNext()) {

				Sequence seq = seqi.nextSequence();
				System.out.println("Sequence length=" + seq.length());
				vbi_id = seq.getName();

				System.out.println("seq name = " + vbi_id);

				SymbolList syms = seq.subList(1, seq.length() - (seq.length())
						% 3);

				Iterator i = syms.iterator();

				while (i.hasNext()) {

					System.out.println("sequence="
							+ ((Symbol) i.next()).toString());
				}

				// System.out.println(syms.getAlphabet());

				// if it is DNA transcribe it to RNA
				if (syms.getAlphabet() == DNATools.getDNA()) {
					syms = DNATools.toRNA(syms);
				}

				if (c == 1) {
					break;
				}
			}
		} catch (Exception ex) {
			ex.printStackTrace();
		}
	}

	public void translate(String fileName, String seqType) {
		String filename = fileName;
		String type = seqType;
		try {
			SymbolTokenization toke = AlphabetManager.alphabetForName(type)
					.getTokenization("token");

			BufferedReader br = new BufferedReader(new FileReader(filename));

			SequenceIterator seqi = RichSequence.IOTools.readFasta(br, toke,
					null);

			// for each sequence
			while (seqi.hasNext()) {

				Sequence seq = seqi.nextSequence();

				// for each frame
				for (int i = 0; i < 3; i++) {
					SymbolList prot;
					Sequence trans;

					// take the reading frame
					// remember that in a SymbolList the first element has
					// index= 1
					// remember that if the length of the list evenly divisible
					// by three an IllegalArgumentException will be thrown

					SymbolList syms = seq.subList(i + 1, seq.length()
							- (seq.length() - i) % 3);

					// if it is DNA transcribe it to RNA
					if (syms.getAlphabet() == DNATools.getDNA()) {
						syms = DNATools.toRNA(syms);
					}

					// output forward translation to STDOUT
					prot = RNATools.translate(syms);

					trans = SequenceTools.createSequence(prot, "", seq
							.getName()
							+ "TranslationFrame: +" + (i + 1),
							Annotation.EMPTY_ANNOTATION);

					// parseAASeq(trans, i+1);

					subStr(trans, i + 1, "forward", seq.getName());

					/*
					 * This method is deprecated since BioJava 1.5
					 * SeqIOTools.writeFasta(System.out, trans);
					 */
					// RichSequence.IOTools.writeFasta(System.out, trans, null);

					// output reverse frame translation to STDOUT
					syms = RNATools.reverseComplement(syms);

					prot = RNATools.translate(syms);
					trans = SequenceTools.createSequence(prot, "", seq
							.getName()
							+ " TranslationFrame: -" + (i + 1),
							Annotation.EMPTY_ANNOTATION);

					subStr(trans, i + 1, "backword", seq.getName());

					/*
					 * This method is deprecated since BioJava 1.5
					 * SeqIOTools.writeFasta(System.out, trans);
					 */
					// RichSequence.IOTools.writeFasta(System.out, trans, null);
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
		}

	}

	public void translate(String fileName) {

		try {
			SymbolTokenization toke = AlphabetManager.alphabetForName("DNA")
					.getTokenization("token");

			BufferedReader br = new BufferedReader(new FileReader(fileName));

			SequenceIterator seqi = RichSequence.IOTools.readFasta(br, toke,
					null);

			String seqName = "";

			// for each sequence
			while (seqi.hasNext()) {

				Sequence seq = seqi.nextSequence();

				SymbolList reverseSyms = seq.subList(1, seq.length());

				reverseSyms = RNATools.reverseComplement(reverseSyms);

				seqName = seq.getName();

				// for each frame
				for (int i = 0; i < 3; i++) {

					forwardTranslate(seq, i);
					backwardTranslate(reverseSyms, i, seqName);

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

	public void forwardTranslate(Sequence seq, int i) throws Exception {
		SymbolList prot;
		Sequence trans;

		SymbolList syms = seq.subList(i + 1, seq.length() - (seq.length() - i)
				% 3);

		// if it is DNA transcribe it to RNA
		if (syms.getAlphabet() == DNATools.getDNA()) {

			syms = DNATools.toRNA(syms);
		}

		// output forward translation to STDOUT
		prot = RNATools.translate(syms);

		trans = SequenceTools.createSequence(prot, "", seq.getName()
				+ "TranslationFrame: +" + (i + 1), Annotation.EMPTY_ANNOTATION);

		subStr(trans, i + 1, "forward", seq.getName());

	}

	public void backwardTranslate(SymbolList reverseSyms, int i, String seqName)
			throws Exception {
		SymbolList prot;
		Sequence trans;

		SymbolList syms = reverseSyms.subList(i + 1, reverseSyms.length()
				- (reverseSyms.length() - i) % 3);

		// if it is DNA transcribe it to RNA
		if (syms.getAlphabet() == DNATools.getDNA()) {
			syms = DNATools.toRNA(syms);
		}

		// output forward translation to STDOUT
		prot = RNATools.translate(syms);

		trans = SequenceTools.createSequence(prot, "", seqName
				+ "TranslationFrame: +" + (i + 1), Annotation.EMPTY_ANNOTATION);

		// parseAASeq(trans, i+1);

		subStr(trans, i + 1, "forward", seqName);

	}

	public static void main(String[] args) {

		Hex hex = new Hex();
		String fileName = "H:/fasta/test.fna";
		String seqType = "DNA";
		hex.translate(fileName, seqType);
		// hex.translate(fileName);
	}

	public void subStr(Sequence seq, int frame, String strand, String seqName) {

		int length = seq.length();
		int start = 1;
		int stop = 1;

		int order = 1;

		System.out.println("seq length= " + length);

		String temp = seq.subStr(1, length);

		int s = 0;
		char t = '|';
		StringBuffer sb = new StringBuffer();

		while (s < length) {
			t = temp.charAt(s);

			if (t == '*') {

				ORF orf = new ORF();

				orf.setStart(start);

				stop = (s) * 3;

				orf.setStop(stop);

				start = stop + 4;

				orf.setFrame(strand + "_" + frame);

				orf.setId(seqName + "_" + orf.getStart() + "_" + orf.getStop()
						+ "_" + orf.getFrame());

				orf.setAaseq(sb.toString());

				orf.setOrder(order);
				orf.setS(s);
				order = order + 1;
				orf.display();
				orfs.add(orf);
				sb = null;
				sb = new StringBuffer();

			} else {
				sb.append(t);
			}
			s = s + 1;
		}

	}

	public Vector<ORF> getOrfs() {
		return orfs;
	}

	public void setOrfs(Vector<ORF> orfs) {
		this.orfs = orfs;
	}

}
