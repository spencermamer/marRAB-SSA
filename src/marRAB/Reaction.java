package marRAB;

import java.util.ArrayList;
import java.util.List;

public class Reaction {

	private List<RxnPair> reactants;
	private List<RxnPair> products;
	
	private double rateConst;
	private double propensity;
	
	public Reaction(double k)
	{
		this.reactants = new ArrayList<RxnPair>();
		this.products = new ArrayList<RxnPair>();
		this.rateConst = k;
	}
	
	public String returnReactants() {
		String s = "";
		for (RxnPair i : reactants) {
			s = s + i.returnMoleculeName() + ":\n\tStoich=" + i.returnStoich() + "\n";
		}
		return s;
	}
	
	public void addReactant(Molecule m, int stoich) {
		reactants.add(new RxnPair(m,stoich));
	}
	
	public void addReactant(Molecule m) {
		addReactant(m,1);
		
	}
	
	public void addProduct(Molecule m, int stoich) {
		products.add(new RxnPair(m, stoich));
	}
	
	public void addProduct(Molecule m) {
		addProduct(m,1);
	}
	
	public double returnPropensity() {
		return propensity;
	}
	
	public double computePropensity() {
		 
		propensity = combinations() * rateConst;
		return propensity;
	}
	
	private double combinations()
	{
		double h = 1.0;
		for(RxnPair reactant : reactants) {
			int m = reactant.returnStoich();
			if(m>1) {
				int num = 1;
				int den = 1;
				for(int times = 0; times<m;times++){
					num *= reactant.returnMolecule().returnNumber() - times;
					den *= m-times;
				}
				h *= num/den;
			} else {
				h *= reactant.returnMolecule().returnNumber();
			}
		}
		
		return h;
	}
	
	public void fireReaction() {
		
		for(RxnPair reactant : reactants) {
			
			if(reactant.returnStoich() > reactant.returnMolecule().returnNumber()) {
				throw new RuntimeException("[ERROR]: Reaction " + this + "failed. Insufficient specie " + reactant.returnMoleculeName()+". Check for improper stoichiometry");
			} else {
				reactant.decrease();
			}
		}
		
		for(RxnPair product : products) {
			product.increase();
		}
	}
}