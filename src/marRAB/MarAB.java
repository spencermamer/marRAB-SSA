package marRAB;


public class MarAB 
{
	public static void main(String[] args)
	{
		double t_end = 1200; // Length of time in minutes
		int imax = 1; // Number of simulations to run
		Model model; // Declare model object
		for(int i = 0; i<imax; i++)	{ // Loop through number of simulations
			model = new Model(t_end, i); // Initialize new model
			model.runModel();
		}
	}
}