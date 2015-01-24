package marRAB;


import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class Model {

	private int logInterval;
	private int simcount;

	private Log lw;

	private double tau;
	private double t;
	private double t_end;
	private Random rand;
	
	private ReactionsList reactionsList;
	private List<Molecule> outputSpecies;
	

	public Model(double t_end, int simcount, int logInterval) {
		this.simcount = (simcount+1);
		this.reactionsList = new ReactionsList();
		this.t_end = t_end;
		this.logInterval = logInterval;
		this.outputSpecies = new ArrayList<Molecule>();
		this.t = 0;
		this.rand = new Random();
		this.lw = new Log(this.outputSpecies, this.simcount);
	}
	
	public Model(double t_end, int simcount) {
		this(t_end, simcount, 5);
	}
	
	public Model(double t_end) { 
		this(t_end, 0);
	}	
	
	public void runModel() {


		// Clear out values that might remain from previous simulations
		outputSpecies.clear();
		reactionsList.empty();
		
		int step = 0; // Initialize loop counter at 0
		t = 0; //Reset t to 0 minutes

		initMarModel_rob(); // Initialize model

		long startTime = System.currentTimeMillis(); //Record time simulation began
		lw.initializeLogFile();
		
		double a0 = reactionsList.sumPropensities(); // Compute initial propensity

		System.out.println("******* Simulation #"+simcount+" beginning. Length = "+t_end+" min");

		while (t < t_end && a0 != 0.0) {
			double r1 = rand.nextDouble();
			double r2 = rand.nextDouble();

			a0 = reactionsList.return_a0();
			tau = ((1.0/a0) * Math.log(1.0/r1));
			t += tau;

			Reaction reaction = reactionsList.determineMu((r2 * a0));
			reaction.fireReaction();

			step++;
			a0 = reactionsList.sumPropensities();
			
			if(step % logInterval == 0) {
				lw.recordState(t);
			}
			
		}
		lw.endLog(t);
		System.out.println("# Execution time: " + ((System.currentTimeMillis() - startTime) / 1000.0) + "s" +
				" Steps: " + step);
	}
	
	private void initMarModel_rob() {

		Molecule a = new Molecule("A", 500); // [Conc]  * (um)^3 /Na = molecules
		Molecule r_2 = new Molecule("R2", 500);
		Molecule p00 = new Molecule ("P00", 1);
		Molecule p01 = new Molecule ("P01", 0);
		Molecule p02 = new Molecule ("P02", 0);
		Molecule p10 = new Molecule ("P10", 0);
		Molecule p11 = new Molecule ("P11", 0);
		Molecule p12 = new Molecule ("P12", 0);
		Molecule m = new Molecule ("M", 0);
		Molecule ruf = new Molecule("Ruf", 0);
		Molecule auf = new Molecule("Auf", 0);
		Molecule r = new Molecule("R", 0);
		Molecule rob = new Molecule ("rob_free", 1000);
		
		Molecule pi0 = new Molecule ("Pi0");
		Molecule pi1 = new Molecule("Pi1");
		Molecule pi2 = new Molecule ("Pi2");
		
		outputSpecies.add(a);
		outputSpecies.add(r_2);
		outputSpecies.add(p00);
		outputSpecies.add(p01);
		outputSpecies.add(p02);
		outputSpecies.add(p10);
		outputSpecies.add(p11);
		outputSpecies.add(p12);
		outputSpecies.add(m);
		outputSpecies.add(ruf);
		outputSpecies.add(auf);
		outputSpecies.add(r);
		outputSpecies.add(rob);
		outputSpecies.add(pi0);
		outputSpecies.add(pi1);
		outputSpecies.add(pi2);
		
		/* Begin initializing reactions
		 * r1 through r14 Promoter dynamics
		 * r15 through r27 Translation/Folding/Transcription
		 * r28 through r32 Degradation
		 */

		Reaction	r1	 = new Reaction(1.8/1500); 
		// a + p00 -> p10
		r1.addReactant(a);
		r1.addReactant(p00);
		r1.addProduct(p10);

		Reaction	r2	 = new Reaction(	1.8	 ); 
		// p10 -> a + p00
		r2.addReactant(p10);
		r2.addProduct(a);
		r2.addProduct(p00);

		Reaction	r3	 = new Reaction(1.8/1500.0/1.5); 
		// a + p01 -> p11
		r3.addReactant(a);
		r3.addReactant(p01);
		r3.addProduct(p11);

		Reaction	r4	 = new Reaction(	1.8	 ); 
		// p11 -> a + p01
		r4.addReactant(p11);
		r4.addProduct(a);
		r4.addProduct(p01);

		Reaction	r5	 = new Reaction(1.8/1500.0/1.5);
		// a + p02 -> p12
		r5.addReactant(p02);
		r5.addReactant(a);
		r5.addProduct(p12);

		Reaction	r6	 = new Reaction(	1.8	 );
		// p12 -> a + p02
		r6.addReactant(p12);
		r6.addProduct(a);
		r6.addProduct(p02);

		Reaction	r7	 = new Reaction(2.0*1.8/150.0);
		// p00 + r2 -> p01
		r7.addReactant(p00);
		r7.addReactant(r_2);
		r7.addProduct(p01);

		Reaction	r8	 = new Reaction(	1.8	 );
		// p01 -> p00 + r2
		r8.addReactant(p01);
		r8.addProduct(p00);
		r8.addProduct(r_2);

		Reaction	r9	 = new Reaction(1.8/150.0);
		// p01 + r2 -> p02
		r9.addReactant(p01);
		r9.addReactant(r_2);
		r9.addProduct(p02);

		Reaction	r10	 = new Reaction(2.0*1.8);
		// p02 -> p01 + r2
		r10.addReactant(p02);
		r10.addProduct(p01);
		r10.addProduct(r_2);

		Reaction r11	 = new Reaction(2.0*1.8/180/1000 );
		// p10 + r2 -> p11
		r11.addReactant(p10);
		r11.addReactant(r_2);
		r11.addProduct(p11);

		Reaction	r12	 = new Reaction(	1.8	 );
		// p11 -> p10 + r2
		r12.addReactant(p11);
		r12.addProduct(p10);
		r12.addProduct(r_2);

		Reaction	r13	 = new Reaction( 1.8/180/1.5);
		// p11 + r2 -> p12
		r13.addReactant(p11);
		r13.addReactant(r_2);
		r13.addProduct(p12);

		Reaction	r14	 = new Reaction(	3.6	 );
		// p12 -> p11 + r2
		r14.addReactant(p12);
		r14.addProduct(p11);
		r14.addProduct(r_2);

		Reaction	r15	 = new Reaction(	0.4	 );
		// p00 -> p00 + m + ruf + auf
		r15.addReactant(p00);
		r15.addProduct(p00);
		r15.addProduct(m);
		r15.addProduct(ruf);
		r15.addProduct(auf);

		Reaction	r16	 = new Reaction(0.4*1.0/800.00 );
		// p01 -> p01 + m + ruf + auf
		r16.addReactant(p01);
		r16.addProduct(p01);
		r16.addProduct(m);
		r16.addProduct(ruf);
		r16.addProduct(auf);

		Reaction	r17	 = new Reaction(	0.4*80.0	 );
		// p10 -> p10 + m + ruf + auf
		r17.addReactant(p10);
		r17.addProduct(p10);
		r17.addProduct(m);
		r17.addProduct(ruf);
		r17.addProduct(auf);

		Reaction	r18	 = new Reaction(	0.4*80.0/800.0	 );
		// p11 -> p11 + m + ruf + auf
		r18.addReactant(p11);
		r18.addProduct(p11);
		r18.addProduct(m);
		r18.addProduct(ruf);
		r18.addProduct(auf);

		Reaction	r19	 = new Reaction(	0.4*80.0/800.0/10.0	 );
		// p12 -> p12 + m + ruf + auf
		r19.addReactant(p12);
		r19.addProduct(p12);
		r19.addProduct(m);
		r19.addProduct(ruf);
		r19.addProduct(auf);		

		Reaction	r20	 = new Reaction(	0.4/800.0/10.0	 );
		// p02 -> p02 + m + ruf + auf
		r20.addReactant(p02);
		r20.addProduct(p02);
		r20.addProduct(m);
		r20.addProduct(ruf);
		r20.addProduct(auf);

		Reaction	r21	 = new Reaction(	0.88	 );
		// m -> m + ruf
		r21.addReactant(m);
		r21.addProduct(m);
		r21.addProduct(r);

		Reaction	r22	 = new Reaction(	6.8	 );
		// m -> m + auf
		r22.addReactant(m);
		r22.addProduct(m);
		r22.addProduct(a);

		Reaction	r23	 = new Reaction(	5	 );
		// auf -> a
		r23.addReactant(auf);
		r23.addProduct(a);

		Reaction	r24	 = new Reaction(	5	 );
		// ruf -> r
		r24.addReactant(ruf);
		r24.addProduct(r);


		Reaction	r25	 = new Reaction(	0.01	 );
		// 2r -> r_2
		r25.addReactant(r,2);
		r25.addProduct(r_2);

		Reaction	r26	 = new Reaction(	0.01/50	 );
		// r_2 -> 2 r
		r26.addReactant(r_2);
		r26.addProduct(r,2);

		Reaction	r27	 = new Reaction(	0.028881133	 );
		// m-> 0
		r27.addReactant(m);

		Reaction	r28	 = new Reaction(	0.028881133	 );
		// auf->0
		r28.addReactant(auf);

		Reaction	r29	 = new Reaction(	0.693147181	 );
		//a->0
		r29.addReactant(a);

		Reaction	r30	 = new Reaction(	0.028881133	 );
		//ruf->0
		r30.addReactant(ruf);

		Reaction	r31	 = new Reaction(	0.028881133	 );
		//r->0
		r31.addReactant(r);

		Reaction	r32	 = new Reaction(	0.028881133	 );
		//r_2->0
		r32.addReactant(r_2);
			
		
		Reaction	r33	 = new Reaction(1.8/1500/2); 
		// rob + p00 -> pi0
		r33.addReactant(rob);
		r33.addReactant(p00);
		r33.addProduct(pi0);

		Reaction	r34	 = new Reaction(	1.8	 ); 
		// p10 -> a + p00
		r34.addReactant(pi0);
		r34.addProduct(rob);
		r34.addProduct(p00);
		
		Reaction 	r35 = new Reaction(1.8/1500.0/1.5/2);
		// rob + p01 -> pi1
		r35.addReactant(rob);
		r35.addReactant(p01);
		r35.addProduct(pi1);
		
		Reaction	r36 = new Reaction(1.8);
		// pi1 -> rob + p01
		r36.addReactant(pi1);
		r36.addProduct(rob);
		r36.addProduct(p01);
		
		Reaction	r37 = new Reaction(1.8/1500.0/1.5/2);
		//  rob + p02 -> pi2
		r37.addReactant(rob);
		r37.addReactant(p02);
		r37.addProduct(pi2);
		
		Reaction 	r38 = new Reaction(1.8);
		// pi2 -> rob + p02
		r38.addReactant(pi2);
		r38.addProduct(rob);
		r38.addProduct(p02);
		
		Reaction	r39 = new Reaction(10*0.4*80.0);
		// pi0 -> pi0 + M + auf + ruf
		r39.addReactant(pi0);
		r39.addProduct(pi0);
		r39.addProduct(m);
		r39.addProduct(auf);
		r39.addProduct(ruf);
		
		Reaction	r40 = new Reaction(10*0.4*80.0/800.0);
		// pi1 -> pi1 + M + auf + ruf
		r40.addReactant(pi1);
		r40.addProduct(pi1);
		r40.addProduct(m);
		r40.addProduct(auf);
		r40.addProduct(ruf);
		
		Reaction	r41 = new Reaction(10*0.4*80.0/800.0/10.0);
		// pi2 -> pi2 + M + auf + ruf
		r41.addReactant(pi2);
		r41.addProduct(pi2);
		r41.addProduct(m);
		r41.addProduct(auf);
		r41.addProduct(ruf);
		
		/*
		 * Insert here reactions where 
		 * 
		 *		1) rob binds the MarA binding spot (P00, P01, P02)-> (Pi0, Pi1, Pi2)
		 *		2) rob unbinds the MarA binding spot 
		 *		3) Transcription of rob-bound sequence
		 *
		 * For now, make the binding 1/2 as strong as for MarA and the transcription twice as fast as for MarA-bound sequence
		 * 
		 * 
		 */


		reactionsList.addReaction(r1);
		reactionsList.addReaction(r2);
		reactionsList.addReaction(r3);
		reactionsList.addReaction(r4);
		reactionsList.addReaction(r5);
		reactionsList.addReaction(r6);
		reactionsList.addReaction(r7);
		reactionsList.addReaction(r8);
		reactionsList.addReaction(r9);
		reactionsList.addReaction(r10);
		reactionsList.addReaction(r11);
		reactionsList.addReaction(r12);
		reactionsList.addReaction(r13);
		reactionsList.addReaction(r14);
		reactionsList.addReaction(r15);
		reactionsList.addReaction(r16);
		reactionsList.addReaction(r17);
		reactionsList.addReaction(r18);
		reactionsList.addReaction(r19);
		reactionsList.addReaction(r20);
		reactionsList.addReaction(r21);
		reactionsList.addReaction(r22);
		reactionsList.addReaction(r23);
		reactionsList.addReaction(r24);
		reactionsList.addReaction(r25);
		reactionsList.addReaction(r26);
		reactionsList.addReaction(r27);
		reactionsList.addReaction(r28);
		reactionsList.addReaction(r29);
		reactionsList.addReaction(r30);
		reactionsList.addReaction(r31);
		reactionsList.addReaction(r32);
		reactionsList.addReaction(r33);
		reactionsList.addReaction(r34);
		reactionsList.addReaction(r35);
		reactionsList.addReaction(r36);
		reactionsList.addReaction(r37);
		reactionsList.addReaction(r38);
		reactionsList.addReaction(r39);
		reactionsList.addReaction(r40);
		reactionsList.addReaction(r41);
		System.out.println("***** marAB initialized");
	}
}



