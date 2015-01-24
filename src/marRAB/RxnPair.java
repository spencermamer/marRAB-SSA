package marRAB;

public class RxnPair {
	private Molecule m;
	private int stoich;
	
	public RxnPair(Molecule m, int stoich) {
		this.m = m;
		this.stoich = stoich;
	}
	
	public RxnPair(Molecule m) {
		this(m, 1);
	}
	
	public Molecule returnMolecule() {
		return m;
	}
	
	public int returnStoich() {
		return stoich;
	}
	
	public String returnMoleculeName() {
		return m.returnName();
	}
	
	public void increase() {
		m.changeNumber(stoich);
	}
	
	public void decrease() {
		m.changeNumber(-stoich);
	}
	
	public void setReagantStoich(int newStoich) {
		stoich = newStoich;
	}
	
}
