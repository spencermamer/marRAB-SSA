package marRAB;

import java.util.List;
import java.util.ArrayList;

public class Compartment {
	
	private String name;
	private double volume;
	private List<Molecule> species;
	
	
	public Compartment(String name) {
		this.name = name;
		this.volume = 1.0;
		
	}
	
	
}
