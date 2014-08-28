package edu.vt.vbi.biojava;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;

public class ORF {
	
	private int start;
	private int stop;
	private String frame;
	private String id;
	private String aaseq;
	private int order;
	
	private int s;
	
	
	public int getStart() {
		return start;
	}
	public void setStart(int start) {
		this.start = start;
	}
	public int getStop() {
		return stop;
	}
	public void setStop(int stop) {
		this.stop = stop;
	}
	public String getFrame() {
		return frame;
	}
	public void setFrame(String frame) {
		this.frame = frame;
	}
	public String getId() {
		return id;
	}
	public void setId(String id) {
		this.id = id;
	}
	public String getAaseq() {
		return aaseq;
	}
	public void setAaseq(String aaseq) {
		this.aaseq = aaseq;
	}
	public int getOrder() {
		return order;
	}
	public void setOrder(int order) {
		this.order = order;
	}
	
	public void display(){
		System.out.println("\n\nid ="+this.id);
		System.out.println("order ="+this.order +" s= "+s);
		System.out.println("seq ="+this.aaseq);
	}
	public int getS() {
		return s;
	}
	public void setS(int s) {
		this.s = s;
	}
	
	public void writeToFile(BufferedWriter w)throws IOException{
		
		w.write(">"+this.id);
		w.write("\n");
		w.write(aaseq);
		w.write("\n");
		w.flush();
		
	}
	
}
