package marRAB;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.List;

public class Log {
	
	private File file;
	private PrintWriter pw;
	List<Molecule> outputSpecies;
	
	// Class constructors
	public Log(List<Molecule> outputSpecies, String fileNameBase, int logCount) {
		this.outputSpecies = outputSpecies;
		this.file = new File(fileNameBase + "_" + logCount + ".txt");
	}
	
	public Log(List<Molecule> outputSpecies, int logCount) {
		this(outputSpecies,("Log_" + System.currentTimeMillis()), logCount);
	}
	
	public Log(List<Molecule> outputSpecies) {
		this(outputSpecies, 0);
	}
	
	
	
	public void initializeLogFile() {
		/* 
		 * Creates PrintWriter instance, creates column titles from Species names, and records initial states
		 * 
		 */
		createPrintWriter();
		createColumnTitles();
		recordState(0.0);
	}
	
	private void createPrintWriter() {
		try {
			pw = new PrintWriter(file);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	private void createColumnTitles() {
		pw.print("Time\t");	
		for(Molecule m : outputSpecies) {
			pw.print(m.returnName() + "\t");
		}
		pw.println();
		pw.flush();
		
	}
	
	public void recordState(double t) {
		pw.print(t+"\t");
		for(Molecule m:outputSpecies) {
			pw.print(m.returnNumber() + "\t");
		}
		pw.println();
		pw.flush();
	}
	
	public void endLog(double t) {
		recordState(t);
		closePrintWriter();
	}
	
	private void closePrintWriter() {
		pw.close();
	}
}
