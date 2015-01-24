package marRAB;

public class Molecule
{
	private String name;
	private int n;
	
	public Molecule(String name) {
		this.name = name;
		this.n = 0;
		}
	
	public Molecule(String name, int n) {
		this.name = name;
		this.n = n;
	}
	

	public String returnName() {
		return name; 
		}
	
	public int returnNumber() {
		return n;
	}
	

	public void changeNumber(int del) {
		
			n += del;
			assert n >= 0;
		
	}
		
	
}