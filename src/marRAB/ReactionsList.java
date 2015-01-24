package marRAB;

import java.util.ArrayList;
import java.util.List;

public class ReactionsList 
{
	private List<Reaction> reactions; // List of the set of reactions describing the reaction
	private double a_0; // Sum of the propensities
	
	public ReactionsList()
	{
		reactions = new ArrayList<Reaction>(); //Initiate
	}
	
	public void init() 
	{
		sumPropensities();
	}
	
	
	public double return_a0()
	{
		return a_0;
	}
	
	
	public void empty()
	{
		this.reactions.clear();
	}
	
	public Reaction addReaction(Reaction r) 
	{
		reactions.add(r);
		return r;
	}
	
	public double sumPropensities()
	{
		a_0 = 0;
		
		//Sums the propensity for each reaction in the reactionsList
		
		for(Reaction r : reactions) { 
			a_0 += r.computePropensity();
		}
		
		return a_0;
	}
	
	public Reaction determineMu(double a) 
	{
		double sum = 0;
		
		for(Reaction r:reactions) 
		{
			sum += r.returnPropensity();
			
			if(sum >= a) 
			{
				
				return r;
			}
		}
		
		throw new RuntimeException("Error: Could find no value of mu; reaction failed for a = " + a);
	}
	
	
}