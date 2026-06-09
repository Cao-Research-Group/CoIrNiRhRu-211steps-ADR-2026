package HEA_fcc_CoIrNiRhRu_211_GitHub;



import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.NumberFormat;

import matsci.Species;
import matsci.engine.monte.metropolis.Metropolis;
import matsci.io.app.log.Status;
import matsci.io.clusterexpansion.PRIM;
import matsci.io.vasp.POSCAR;
import matsci.location.symmetry.operations.SpaceGroup;
import matsci.location.symmetry.operations.SymmetryOperation;
import matsci.structure.Structure;
import matsci.structure.decorate.DecorationTemplate;
import matsci.structure.decorate.function.ce.ClusterExpansion;
import matsci.structure.decorate.function.ce.FastAppliedCE;
import matsci.structure.decorate.function.ce.clusters.ClusterGroup;
import matsci.structure.superstructure.SuperStructure;

public class runMonteCanonical {

	public static String ROOT_DIR = "ce/CoIrNiRhRu-211/";
    public static String GROUND_DIR=ROOT_DIR + "/groundstates/";
    public static String PRIM_FILE = ROOT_DIR + "/PRIM.vasp";
    public static String CLUSTER_DIR = ROOT_DIR + "/clusters/";



    public static void main(String[] args) {

        Status.setLogLevelDetail();

        ClusterExpansion ce = getPreFittedCE();

    	String[] abc = new String[args.length];
        for(int i = 0; i < args.length & i < 100; i++){
               abc[i]=args[i];
        }

        runMonteCanonical_Combined(ce, 5, 0, 0, 16, "2000K-equimolar-numpass=5000", "Sym_equimolar-prim=80.vasp", 5000, 2000, 980);
        System.out.println("OK");

    }

	  /**
	   *
	   *  runMonte; GC, grand canonical ensemble; chemPot, set chemical potential
	   * @param ce Cluter expansions
	   * @param chemPotPtNi, the chemical potentail difference between Pt and Ni
	   * @param chemPotOxygen, the chemical potential difference between oxygen and vacancy
	   * @param size1, the x expansion of supercell
	   * @param size2, the y expansion of supercell
	   * @return
	   */
	  public static int[] runMonteCanonical_Combined(ClusterExpansion ce, int sizex1, int sizex2, int sizey1, int sizey2, String fileDir, String initialStr, int numPass,double startTemp2,int endTemp2) {

	      NumberFormat numberFormat2 = NumberFormat.getNumberInstance();
	      numberFormat2.setMaximumFractionDigits(4);

	      Status.basic("\n\nComposition profiles for CoIrNiRhRu under 3000K to 1000K.");
	      

	      String fileDirTemp = GROUND_DIR + fileDir + "/" + "intermediateStrs/";

	      File dirFileTemp  = new  File(fileDirTemp);
	      if ( ! (dirFileTemp.exists())  &&   ! (dirFileTemp.isDirectory())) {
	              dirFileTemp.mkdirs();
	      }

	      double startTemp = startTemp2;
	      double endTemp = endTemp2;
	      double tempIncrement = 20;
	      double concentrationCo = 0;
	      double concentrationIr = 0;
	      double concentrationNi = 0;
	      double concentrationRh = 0;
	      double concentrationRu = 0;
	      double concentrationN = 0;
	      //double concentrationVac = 0;
	      double calculatedEnergy = 0;
	      int numPrims = 0;
	      int numPasses = numPass;
	      System.out.println("the number of numPasses : " + numPasses);

	      /*
	       * chemical potential of all species
	       */

	      double chemPotCo = 0.0;
	      double chemPotIr = 0.0;
	      double chemPotNi = 0.0;
	      double chemPotRh = 0.0;
	      double chemPotRu = 0.0;
	      double chemPotN = 0.0;

	      Status.detail(ce.getClusterGroup(2).getECI(0) + " ");
	      Status.detail(ce.getClusterGroup(3).getECI(0) + " ");
	      Status.detail(ce.getClusterGroup(4).getECI(0) + " ");


	      int[][] superToDirect = new int[][] {
	            {sizex1, sizex2},
	            {sizey1, sizey2}
	      };

	      double kb = 8.6173423E-5;

	      FastAppliedCE appliedCE = new FastAppliedCE(ce, superToDirect);
	      appliedCE.activateAllGroups();

	      System.out.println("the number of sites is: " + appliedCE.numSigmaSites());

	      /*
	       * set chemical potential for different species, can play around to get how it will effect the "ground" state.
	       * pay attention to the position where the chemical potential should be set.
	       */

	      appliedCE.setChemPot(Species.cobalt, chemPotCo);
	      appliedCE.setChemPot(Species.iridium, chemPotIr);
	      appliedCE.setChemPot(Species.nickel, chemPotNi);
	      appliedCE.setChemPot(Species.rhodium, chemPotRh);
	      appliedCE.setChemPot(Species.ruthenium, chemPotRu);
	      appliedCE.setChemPot(Species.nitrogen, chemPotN);

	      Metropolis metropolis = new Metropolis();
	      metropolis.setNumIterations(1000);
	      appliedCE.decorateWithSpecies(Species.ruthenium);


	      System.out.println("print out the initial structure before decoration");
	      Structure structure12 = appliedCE.getStructure();
	      POSCAR outfile= new POSCAR(structure12);
	      outfile.writeFile(fileDirTemp + "groundState_1.vasp");
	      outfile.writeVICSOutFile(fileDirTemp + "groundState_1.out");
	      POSCAR sourcePoscar = new POSCAR(GROUND_DIR + initialStr);
	      sourcePoscar.setVectorPeriodicity(2, false);
	      Structure sourceStr = new Structure(sourcePoscar);


	      SpaceGroup space1=sourceStr.getDefiningSpaceGroup();

	      //get the symmetrical operators corresponding to space group
	      SymmetryOperation[] operation1=space1.getOperations();

	      SymmetryOperation mirrorOp = operation1[0];
	      System.out.println("the number of symmetrical operations is : " + operation1.length);

	      for(int i=0; i< operation1.length; i++){
	          SymmetryOperation ope1 = operation1[i];
	          double[][] abc =ope1.getOriginCenteredPointOperator(ope1.getBasis());
	          if( (abc[2][2] + 1) < 0.001 && (abc[2][2] + 1) > -0.001 ) {
	              System.out.println("operation  " + i + " is the mirror operations");
	              mirrorOp = operation1[i];
	              break;
	          }
	      }
	


	      appliedCE.decorateFromStructure(sourceStr);

	      System.out.println("the number of sites is: " + appliedCE.numSigmaSites());

	      System.out.println("the initial structure has decorated from a source structure, speed the calculation.");

	      int numCo = sourceStr.numDefiningSitesWithSpecies(Species.cobalt);
	      int numIr = sourceStr.numDefiningSitesWithSpecies(Species.iridium);
	      int numNi = sourceStr.numDefiningSitesWithSpecies(Species.nickel);
	      int numRh = sourceStr.numDefiningSitesWithSpecies(Species.rhodium);
	      int numRu = sourceStr.numDefiningSitesWithSpecies(Species.ruthenium);
	      int numN = sourceStr.numDefiningSitesWithSpecies(Species.nitrogen);
	      int numVac = appliedCE.numSigmaSites() - numCo - numIr - numNi - numRh - numRu - numN;

	      concentrationCo = (double)(numCo) / (numCo + numIr + numNi + numRh + numRu);
	      concentrationIr = (double)(numIr) / (numCo + numIr + numNi + numRh + numRu);
	      concentrationNi = (double)(numNi) / (numCo + numIr + numNi + numRh + numRu);
	      concentrationRh = (double)(numRh) / (numCo + numIr + numNi + numRh + numRu);
	      concentrationRu = (double)(numRu) / (numCo + numIr + numNi + numRh + numRu);
	      concentrationN = (double)(numN) / (numN + numVac);
	      //print out the initial structure before metropolis algorithm
	      SuperStructure groundState2 = appliedCE.getSuperStructure();

	      POSCAR outfile0= new POSCAR(groundState2);
	      outfile0.writeFile(fileDirTemp + "groundState_0.vasp");
	      outfile0.writeVICSOutFile(fileDirTemp + "groundState_0.out");


	      if ((endTemp - startTemp) / tempIncrement < 0) {
	         tempIncrement *= -1;
	      }
	      int numSteps = (int) Math.floor((endTemp - startTemp) / tempIncrement);
	      if(startTemp==Double.MAX_VALUE) numSteps=1;

	      try {

	          FileWriter fw = new FileWriter(fileDir + "Temp-AdsE-TOF.txt" , false);
	          BufferedWriter bw = new BufferedWriter(fw);
	          bw.write("Temp	NumCo	NumIr	NumNi	NumRh	NumRu	NumN	formE	formCV	avgAdsE	avgTOF	avgAdsE2	avgTOF2\r\n");

	          FileWriter fw1 = new FileWriter(fileDir + "SublatticeProfile.txt", false);
	          BufferedWriter bw1 = new BufferedWriter(fw1);
	          
	          FileWriter fw2 = new FileWriter(fileDir + "AverageAdsEnergy.txt", false);
	          BufferedWriter bw2 = new BufferedWriter(fw2);
	          
	          FileWriter fw3 = new FileWriter(fileDir + "AverageTOF.txt", false);
	          BufferedWriter bw3 = new BufferedWriter(fw3);
	          
	          FileWriter fw4 = new FileWriter(fileDir + "SublatticeProfile-1st.txt", false);
	          BufferedWriter bw4 = new BufferedWriter(fw4);

	          FileWriter fw5 = new FileWriter(fileDir + "SublatticeProfile-2nd.txt", false);
	          BufferedWriter bw5 = new BufferedWriter(fw5);

	          FileWriter fw6 = new FileWriter(fileDir + "SublatticeProfile-3rd.txt", false);
	          BufferedWriter bw6 = new BufferedWriter(fw6);

	          FileWriter fw7 = new FileWriter(fileDir + "SublatticeProfile-4th.txt", false);
	          BufferedWriter bw7 = new BufferedWriter(fw7);

	          FileWriter fw8 = new FileWriter(fileDir + "SublatticeProfile-5th.txt", false);
	          BufferedWriter bw8 = new BufferedWriter(fw8);
	          
	          FileWriter fw9 = new FileWriter(fileDir + "SublatticeProfile-6th.txt", false);
	          BufferedWriter bw9 = new BufferedWriter(fw9);
	          
	          FileWriter fw10 = new FileWriter(fileDir + "SublatticeProfile-7th.txt", false);
	          BufferedWriter bw10 = new BufferedWriter(fw10);
	          
	          FileWriter fw11 = new FileWriter(fileDir + "SublatticeProfile-8th.txt", false);
	          BufferedWriter bw11 = new BufferedWriter(fw11);
	          
	          FileWriter fw12 = new FileWriter(fileDir + "SublatticeProfile-9th.txt", false);
	          BufferedWriter bw12 = new BufferedWriter(fw12);
	          
	          FileWriter fw13 = new FileWriter(fileDir + "SublatticeProfile-10th.txt", false);
	          BufferedWriter bw13 = new BufferedWriter(fw13);

	          FileWriter fw14 = new FileWriter(fileDir + "SublatticeProfile-11th.txt", false);
	          BufferedWriter bw14 = new BufferedWriter(fw14);
	          
	          FileWriter fw15 = new FileWriter(fileDir + "SublatticeProfile-Adsorption.txt", false);
	          BufferedWriter bw15 = new BufferedWriter(fw15);

	          for (int stepNum = 0; stepNum < numSteps; stepNum++) {

	              double temp = startTemp + stepNum * tempIncrement;
	

	              metropolis.setTemp(kb * temp);
	              System.out.println();
	              /*
	               * mix ensemble ensemble, the middle layer is fixed.
	               */
	              DecorationMSymCanonicalMidFreezeManager manager = new DecorationMSymCanonicalMidFreezeManager(appliedCE, mirrorOp);
	              System.out.println("the manager used here is: DecorationMSymCanonicalMidFreezeManager, mirror sym and freeze middle layer Canonical");

	              CombinedRecorder recorder = new CombinedRecorder(manager, appliedCE);

	              int numPassTemp = numPasses;


	              // just for symmetrical Monte Carlo simulations
	              // Equilibration
	              metropolis.setNumIterations(numPassTemp * appliedCE.numSigmaSites()/10*9);
	              metropolis.runBasic(manager);

	              // Run for real
	              metropolis.setNumIterations(numPassTemp * appliedCE.numSigmaSites()/10);
	              metropolis.run(manager, recorder);

	              System.out.println("Temperature: " + metropolis.getTemp());
	              System.out.println("Trigger ratio: " + recorder.getTriggerRatio());
	              System.out.println("Average Energy: " + recorder.getAverageValue());
	              System.out.println("Average Adsorption Energy: " + recorder.getAverageAdsEnergy());
	              System.out.println("Average TOF: " + recorder.getAverageActivity());

	              /*
	               * print out the ground state for every temperature step
	               *
	               */
	              SuperStructure groundState3 = appliedCE.getSuperStructure();
	              Structure primGroundState3=groundState3.findPrimStructure().getCompactStructure();
	              int numCoTem = primGroundState3.numDefiningSitesWithSpecies(Species.cobalt);
	              int numIrTem = primGroundState3.numDefiningSitesWithSpecies(Species.iridium);
	              int numNiTem = primGroundState3.numDefiningSitesWithSpecies(Species.nickel);
	              int numRhTem = primGroundState3.numDefiningSitesWithSpecies(Species.rhodium);
	              int numRuTem = primGroundState3.numDefiningSitesWithSpecies(Species.ruthenium);
	              int numNTem = primGroundState3.numDefiningSitesWithSpecies(Species.nitrogen);
	              POSCAR outfile4= new POSCAR(primGroundState3);
	              outfile4.writeFile(fileDirTemp + "groundState." + stepNum + ".vasp");
	              outfile4.writeVICSOutFile(fileDirTemp + "groundState." + stepNum  + ".out");
	              

	             calculatedEnergy = recorder.getAverageValue();
	             double calculatedEnergy2 = recorder.getAverageSquaredValue();
		         double avgAdsE = recorder.getAverageAdsEnergy();
	             double avgActi = recorder.getAverageActivity();
	             double avgAdsE2 = recorder.getAverageAdsEnergySq();
	             double avgActi2 = recorder.getAverageActivitySq();
	             double FormESq = calculatedEnergy2 - calculatedEnergy * calculatedEnergy;
	             double FormECV = FormESq / (kb * temp * temp);
	             double AdsESq = avgAdsE2 - avgAdsE * avgAdsE;
	             double AdsEcV = AdsESq / (kb * temp * temp);
	             AdsEcV /= appliedCE.numSigmaSites();
	             double ActiSq = avgActi2 - avgActi * avgActi;
	             double ActicV = ActiSq / (kb * temp * temp);
	             ActicV /= appliedCE.numSigmaSites();
	             System.out.println("Temp: " + temp + "\tNumPass:" + numPassTemp + "\tConcentration: " + concentrationCo + " " + concentrationIr + " " + concentrationNi + " " + concentrationRh + " " + concentrationRu + " " + concentrationN  + "\tformE: " + calculatedEnergy + "\tAdsE: " + avgAdsE + "\tcv: " + AdsEcV + "\tTOF: " + avgActi + "\tcv: " + ActicV);
	             System.out.println("the num of steps is: " + stepNum);

	             String stringNum = temp + "     " + numCoTem + "     " + numIrTem + "     " + numNiTem + "     " + numRhTem + "     " + numRuTem + "     " + numNTem + "     " + calculatedEnergy +  "     " + FormECV + "     " + avgAdsE + "     " + avgActi + "     " + AdsEcV + "     " + ActicV + "\r\n";
	             bw.write(stringNum);
	             
	             bw1.write("the temperature is: " + temp +"\n");
	             double[][] sublatticeCount;
	             DecorationTemplate decorTemplate = appliedCE.getHamiltonian();
	             sublatticeCount = new double[decorTemplate.numSublattices()][];

	             for(int sublatticeIndex = 0; sublatticeIndex < decorTemplate.numSublattices(); sublatticeIndex++ ){
	                 int numStates = decorTemplate.getAllowedSpeciesForSublattice(sublatticeIndex).length;
	                 sublatticeCount[sublatticeIndex] = new double[numStates];
	                 int[] sigmaSiteForSublattice = appliedCE.getSigmaIndicesForSublattice(sublatticeIndex);
	                 int numSitesForSublattice = sigmaSiteForSublattice.length;

	                 bw1.write("the Lattice index of  " + sublatticeIndex + "     ");
	                 bw1.write("numOfSites=" + numSitesForSublattice + "     ");
	                 for(int state = 0; state < numStates; state++){
	                      for (int siteNum =0; siteNum < numSitesForSublattice; siteNum++){
	                          sublatticeCount[sublatticeIndex][state] += recorder.getAverageStateCount(sigmaSiteForSublattice[siteNum],state);
	                      }
	                      sublatticeCount[sublatticeIndex][state] = sublatticeCount[sublatticeIndex][state]/numSitesForSublattice;
	                      bw1.write(sublatticeCount[sublatticeIndex][state] + "     ");
	                 }
	                 bw1.write("\n");
	                 bw1.write("\r\n");
	                 
	                 switch (sublatticeIndex) {
                      case 0:{
                          bw4.write("temp " + temp + "     ");
                          bw4.write("numOfSites=" + numSitesForSublattice + "     "); //sitesLattice.length=5
                          bw4.write(sublatticeCount[sublatticeIndex][0] + "     ");
                          bw4.write(sublatticeCount[sublatticeIndex][1] + "     ");
                          bw4.write(sublatticeCount[sublatticeIndex][2] + "     ");
                          bw4.write(sublatticeCount[sublatticeIndex][3] + "     ");
                          bw4.write(sublatticeCount[sublatticeIndex][4] + "     ");
                          bw4.write("\r\n");
                          break;
                      }
                      case 1:{
                          bw5.write("temp " + temp + "     ");
                          bw5.write("numOfSites=" + numSitesForSublattice + "     "); //sitesLattice.length=5
                          bw5.write(sublatticeCount[sublatticeIndex][0] + "     ");
                          bw5.write(sublatticeCount[sublatticeIndex][1] + "     ");
                          bw5.write(sublatticeCount[sublatticeIndex][2] + "     ");
                          bw5.write(sublatticeCount[sublatticeIndex][3] + "     ");
                          bw5.write(sublatticeCount[sublatticeIndex][4] + "     ");
                          bw5.write("\r\n");
                          break;
                      }
                      case 2:{
                          bw6.write("temp " + temp + "     ");
                          bw6.write("numOfSites=" + numSitesForSublattice + "     "); //sitesLattice.length=5
                          bw6.write(sublatticeCount[sublatticeIndex][0] + "     ");
                          bw6.write(sublatticeCount[sublatticeIndex][1] + "     ");
                          bw6.write(sublatticeCount[sublatticeIndex][2] + "     ");
                          bw6.write(sublatticeCount[sublatticeIndex][3] + "     ");
                          bw6.write(sublatticeCount[sublatticeIndex][4] + "     ");
                          bw6.write("\n");
                          break;
                      }
                      case 3:{
                          bw7.write("temp " + temp + "     ");
                          bw7.write("numOfSites=" + numSitesForSublattice + "     "); //sitesLattice.length=5
                          bw7.write(sublatticeCount[sublatticeIndex][0] + "     ");
                          bw7.write(sublatticeCount[sublatticeIndex][1] + "     ");
                          bw7.write(sublatticeCount[sublatticeIndex][2] + "     ");
                          bw7.write(sublatticeCount[sublatticeIndex][3] + "     ");
                          bw7.write(sublatticeCount[sublatticeIndex][4] + "     ");
                          bw7.write("\r\n");
                          break;
                      }
                      case 4:{
                          bw8.write("temp " + temp + "     ");
                          bw8.write("numOfSites=" + numSitesForSublattice + "     "); //sitesLattice.length=5
                          bw8.write(sublatticeCount[sublatticeIndex][0] + "     ");
                          bw8.write(sublatticeCount[sublatticeIndex][1] + "     ");
                          bw8.write(sublatticeCount[sublatticeIndex][2] + "     ");
                          bw8.write(sublatticeCount[sublatticeIndex][3] + "     ");
                          bw8.write(sublatticeCount[sublatticeIndex][4] + "     ");
                          bw8.write("\r\n");
                          break;
                      }
                      case 5:{
                          bw9.write("temp " + temp + "     ");
                          bw9.write("numOfSites=" + numSitesForSublattice + "     "); //sitesLattice.length=5
                          bw9.write(sublatticeCount[sublatticeIndex][0] + "     ");
                          bw9.write(sublatticeCount[sublatticeIndex][1] + "     ");
                          bw9.write(sublatticeCount[sublatticeIndex][2] + "     ");
                          bw9.write(sublatticeCount[sublatticeIndex][3] + "     ");
                          bw9.write(sublatticeCount[sublatticeIndex][4] + "     ");
                          bw9.write("\r\n");
                          break;
                      }
                      case 6:{
                          bw10.write("temp " + temp + "     ");
                          bw10.write("numOfSites=" + numSitesForSublattice + "     "); //sitesLattice.length=5
                          bw10.write(sublatticeCount[sublatticeIndex][0] + "     ");
                          bw10.write(sublatticeCount[sublatticeIndex][1] + "     ");
                          bw10.write(sublatticeCount[sublatticeIndex][2] + "     ");
                          bw10.write(sublatticeCount[sublatticeIndex][3] + "     ");
                          bw10.write(sublatticeCount[sublatticeIndex][4] + "     ");
                          bw10.write("\r\n");
                          break;
                      }
                      case 7:{
                          bw11.write("temp " + temp + "     ");
                          bw11.write("numOfSites=" + numSitesForSublattice + "     "); //sitesLattice.length=5
                          bw11.write(sublatticeCount[sublatticeIndex][0] + "     ");
                          bw11.write(sublatticeCount[sublatticeIndex][1] + "     ");
                          bw11.write(sublatticeCount[sublatticeIndex][2] + "     ");
                          bw11.write(sublatticeCount[sublatticeIndex][3] + "     ");
                          bw11.write(sublatticeCount[sublatticeIndex][4] + "     ");
                          bw11.write("\r\n");
                          break;
                      }
                      case 8:{
                          bw12.write("temp " + temp + "     ");
                          bw12.write("numOfSites=" + numSitesForSublattice + "     "); //sitesLattice.length=5
                          bw12.write(sublatticeCount[sublatticeIndex][0] + "     ");
                          bw12.write(sublatticeCount[sublatticeIndex][1] + "     ");
                          bw12.write(sublatticeCount[sublatticeIndex][2] + "     ");
                          bw12.write(sublatticeCount[sublatticeIndex][3] + "     ");
                          bw12.write(sublatticeCount[sublatticeIndex][4] + "     ");
                          bw12.write("\r\n");
                          break;
                      }
                      case 9:{
                          bw13.write("temp " + temp + "     ");
                          bw13.write("numOfSites=" + numSitesForSublattice + "     "); //sitesLattice.length=5
                          bw13.write(sublatticeCount[sublatticeIndex][0] + "     ");
                          bw13.write(sublatticeCount[sublatticeIndex][1] + "     ");
                          bw13.write(sublatticeCount[sublatticeIndex][2] + "     ");
                          bw13.write(sublatticeCount[sublatticeIndex][3] + "     ");
                          bw13.write(sublatticeCount[sublatticeIndex][4] + "     ");
                          bw13.write("\r\n");
                          break;
                      }
                      case 10:{
                    	  bw14.write("temp " + temp + "     ");
                          bw14.write("numOfSites=" + numSitesForSublattice + "     "); //sitesLattice.length=5
                          bw14.write(sublatticeCount[sublatticeIndex][0] + "     ");
                          bw14.write(sublatticeCount[sublatticeIndex][1] + "     ");
                          bw14.write(sublatticeCount[sublatticeIndex][2] + "     ");
                          bw14.write(sublatticeCount[sublatticeIndex][3] + "     ");
                          bw14.write(sublatticeCount[sublatticeIndex][4] + "     ");
                          bw14.write("\r\n");
                          break;
                      }
                      case 11:{
                          // nitrogen binding sites
                          //break;
                    	  bw15.write("temp " + temp + "     ");
                          bw15.write("numOfSites=" + numSitesForSublattice + "     "); //sitesLattice.length=2
                          bw15.write(sublatticeCount[sublatticeIndex][0] + "     ");
                          bw15.write(sublatticeCount[sublatticeIndex][1] + "     ");
                          bw15.write("\r\n");
                          break;
                      }
                      default:
	                 }
	              }
	             
	              bw2.write("the temperature is: " + temp +"\n");
	              System.out.println("the activity distribution is: ");
	              bw2.write("the adsorption energy distribution is: " + "\n");
	              for (int varNum = 0; varNum < recorder.getRecorderEads40Bars_AdsE().length; varNum++) {
	                    System.out.print(recorder.getRecorderEads40Bars_AdsE()[varNum] + "	");
			            bw2.write(recorder.getRecorderEads40Bars_AdsE()[varNum]+"	");
	              }
	              
	              bw2.write("\n");
	              bw2.write("the average adsorption energy is: " + recorder.getAverageAdsEnergy() + "\n");
                  bw2.write("\n");
                  bw2.write("\r\n");
                  
                  
	              bw3.write("the temperature is: " + temp +"\n");
	              System.out.println("\nthe TOF distribution is: ");
	              bw3.write("the TOF distribution is: " + "\n");
	              for (int varNum = 0; varNum < recorder.getRecorderEads40Bars_Acti().length; varNum++) {
	                    System.out.print(recorder.getRecorderEads40Bars_Acti()[varNum] + "	");
			            bw3.write(recorder.getRecorderEads40Bars_Acti()[varNum]+"	");
	              }
	              
	              bw3.write("\n");
	              bw3.write("the average TOF is: " + recorder.getAverageActivity() + "\n");
                  bw3.write("\n");
                  bw3.write("\r\n");
	          }

	          bw.flush();
	          bw.close();
	          fw.close();

	          bw1.flush();
	          bw1.close();
	          fw1.close();

	          bw2.flush();
	          bw2.close();
	          fw2.close();
	          
	          bw3.flush();
	          bw3.close();
	          fw3.close();

	          bw4.flush();
	          bw4.close();
	          fw4.close();

	          bw5.flush();
	          bw5.close();
	          fw5.close();

	          bw6.flush();
	          bw6.close();
	          fw6.close();

	          bw7.flush();
	          bw7.close();
	          fw7.close();

	          bw8.flush();
	          bw8.close();
	          fw8.close();

	          bw9.flush();
	          bw9.close();
	          fw9.close();

	          bw10.flush();
	          bw10.close();
	          fw10.close();

	          bw11.flush();
	          bw11.close();
	          fw11.close();

	          bw12.flush();
	          bw12.close();
	          fw12.close();

	          bw13.flush();
	          bw13.close();
	          fw13.close();

	          bw14.flush();
	          bw14.close();
	          fw14.close();

	          bw15.flush();
	          bw15.close();
	          fw15.close();
	      }

	      catch (IOException e) {
	          e.printStackTrace();
	      }

	      Structure groundState = appliedCE.getSuperStructure().findPrimStructure().getCompactStructure();
	      int numCoF = groundState.numDefiningSitesWithSpecies(Species.cobalt);
	      int numIrF = groundState.numDefiningSitesWithSpecies(Species.iridium);
	      int numNiF = groundState.numDefiningSitesWithSpecies(Species.nickel);
	      int numRhF = groundState.numDefiningSitesWithSpecies(Species.rhodium);
	      int numRuF = groundState.numDefiningSitesWithSpecies(Species.ruthenium);
	      int numNF = groundState.numDefiningSitesWithSpecies(Species.nitrogen);
	      
	      numPrims = appliedCE.numPrimCells();
	      POSCAR outfile2= new POSCAR(groundState);
	      outfile2.setDescription("calculated avgE: " + calculatedEnergy + "chemPot of CoIrNiRhRu: " + numberFormat2.format(chemPotCo) + "," + numberFormat2.format(chemPotIr) + "," + numberFormat2.format(chemPotNi) + "," + numberFormat2.format(chemPotRh) + "," + numberFormat2.format(chemPotRu));
	      outfile2.writeFile(fileDir + "prim=" + numPrims + "-" + sizex1 + "_" + sizex2 + "_" + sizey1 + "_" + sizey2 + "."+"groundState." + numCoF + "." + numIrF + "." + numNiF + "." + numRhF + "." + numRuF + "-" + numNF + ".vasp");
	      outfile2.writeVICSOutFile(fileDir + "prim=" + numPrims + "-" + sizex1 + "_" + sizex2 + "_" + sizey1 + "_" + sizey2 + "."+"groundState." + numCoF + "." + numIrF + "." + numNiF + "." + numRhF + "." + numRuF + "-" + numNF + ".out");

	      String fileDir3 = GROUND_DIR + fileDir + "/";
	      outfile2.writeFile(fileDir3 + numberFormat2.format(chemPotCo) + "." + numberFormat2.format(chemPotIr) + "." + numberFormat2.format(chemPotNi) + "." + numberFormat2.format(chemPotRh) + "." + numberFormat2.format(chemPotRu) + "-" + sizex1 + "_" + sizex2 + "_" + sizey1 + "_" + sizey2 + "."+"groundState." +numCoF + "." + numIrF + "." + numNiF + "." + numRhF + "." + numRuF + "-" + numNF + ".vasp");
	      outfile2.writeVICSOutFile(fileDir3 + numberFormat2.format(chemPotCo) + "." + numberFormat2.format(chemPotIr) + "." + numberFormat2.format(chemPotNi) + "." + numberFormat2.format(chemPotRh) + "." + numberFormat2.format(chemPotRu) + "-" + sizex1 + "_" + sizex2 + "_" + sizey1 + "_" + sizey2 + "."+"groundState." +numCoF + "." + numIrF + "." + numNiF + "." + numRhF + "." + numRuF + "-" + numNF + ".out");
	      

	      int[] numCoIrNiRhRu =  new int[6];
	       numCoIrNiRhRu[0] = numCoF;
	       numCoIrNiRhRu[1] = numIrF;
	       numCoIrNiRhRu[2] = numNiF;
	       numCoIrNiRhRu[3] = numRhF;
	       numCoIrNiRhRu[4] = numRuF;
	       numCoIrNiRhRu[5] = numNF;

	       return numCoIrNiRhRu;
	  } // end of runMonte method

	   
	  public static ClusterExpansion getPreFittedCE(){		//Stores the ECIs of fitted CE

      ClusterExpansion ce = buildCE(); 
      ClusterGroup[] activeGroups = ce.getAllGroups(true);

      double[] eci = new double[]{
    		  5.147579154817643, 0.013619223496654991, 0.06483149223248909, -0.07105080558508876, -0.02500006539157595, -0.006148169765344559, 0.030180456383057063, -0.04031660228726169, -0.021572996738573647, -0.006609674048189693, 0.061554680533084466, 0.004619850377747151, 0.05530158506709745, -0.020538044538672702, 0.011480886795605284, -0.0039050429188954155, -0.019985299710650772, -0.007960083874459192, 0.0020228066742498246, 0.0033483589638845894, -0.016985711161448182, 0.005790086347704635, 0.011662714259537329, 0.006986043164220189, 0.022544523844196078, 0.009926356121532934, 0.026841911582628435, 0.019426427480411047, 0.029864708185306753, 0.005839026220529956, 0.020113664940191792, 0.0277977021132173, 0.0304279296748576, 0.003075547782145254, 0.021626709894212504, 0.037319954888248186, 0.03573568592423107, 0.005538013178460619, 0.021625980061666873, 0.020493109467326107, 0.03251399236112079, -0.013623293451343992, 0.020739691765631577, 0.025942490617950204, 0.03380095667408338, -0.11946580126103414, -1.1175555212384768E-4, 6.118378219674358E-4, -0.0014426231469845598, -5.931066020941811E-4, 9.04129775355754E-4, 8.180447861976514E-4, 0.0019912934729929303, 0.002792947046548735, -0.0017260465066296104, 0.0017896010289272012, -3.6541950472749034E-4, 0.0015879005070892913, -4.394422283754114E-4, 8.30987876790527E-4, 0.0016099152961712299, 0.005202876324977954, -0.001194151963351305, 2.8415073029953953E-4, -9.175047817636204E-5, 1.1555352456149603E-4, 2.4976756576372723E-4, 8.122551688447826E-4, 7.423561558846624E-6, 0.0015046925392092135, -3.471346691542661E-4, -3.919366151689794E-4, -5.014989250579534E-4, 1.369971856812973E-4, -2.5689122844993007E-4, 8.420106429262494E-4, 1.813329457249881E-5, 0.0029229925766729864, 0.021899266700758716, -0.03439690549850845, 0.015414915497271135, -0.011623775378758784, -1.3682550790088306E-4, 0.0012648505580906654, -0.002893066879602951, 0.002682576031505051, 0.001415216946745324, 0.0017588303130646158, 0.0015659432129107988, -0.0025753054575208835, 0.001370139256695387, 0.0033146266106999916, -0.0012774061999610122, -0.0013974326920691798, -0.0029640253132107203, -1.2602199965661837E-4, 0.0022160603345690935, 1.5910911467823235E-4, 6.439125538066787E-4, 0.0018211351833197888, -0.004066975223735199, 0.0013844178829295877, -6.091086308290744E-4, 0.001841144863240783, -3.1474766451390363E-4, -3.626332618552805E-5, 0.0011267838663051037, 0.0035217365978960392, -0.0021091523062250353, 0.001245980818304009, -0.0013712671254124753, 0.0012034274797352042, 0.0014839234610473486, 0.001852644249671236, 4.689991509537685E-4, -1.8466488047460768E-4, -1.8611118755911286E-4, 0.0011446048291283576, 8.098855782621623E-4, 4.3610676863263705E-5, -0.0010898978567150358, -2.962699340063223E-4, 9.095950616810292E-4, 4.6404850882667947E-4, -0.0010053231544546442, 2.177432653998712E-4, -0.002932658943422406, -9.128506547461462E-4, 0.001103622662820104, 2.0597810961222674E-4, 0.002451813528819965, 0.002613629040890017, -0.003036393670596933, 0.003202130356688564, -0.00203271462651883, -4.315514704868683E-5, -8.2428838546031E-4, 0.0030076880139824113, 0.0019111305591243254, 0.0022372160131235285, 3.6954135333632475E-4, 4.173506405594988E-4, -0.005032670517777218, 6.896071329103904E-5, 9.093306788976515E-4, 0.0021773533543489937, 8.785422704902515E-4, 0.002722151907600797, -0.003837979778094985, 0.0015358571451040098, -0.002084249246304702, 3.851229858007855E-4, 6.736506107517886E-4, 0.003052540512029439, -2.692759358822825E-4, 0.004951975077555841, 2.672050692499586E-4, 3.0124019093131046E-4, -5.188149104082642E-4, 6.236360497585546E-4, 6.761357905305811E-4, 7.024411750418838E-4, 8.156779062679477E-5, -5.469677427869248E-4, -0.0011707793331063787, -4.436775832525157E-5, 0.001942940179272591, 0.0011834967485915416, -7.012612648913861E-5, 6.434716258280015E-4, 0.0021482766492958837, 0.0028383085532652065, -1.0912367052162003E-4, 6.494169846817687E-4, -0.0017210215116721778, 6.03501079769824E-4, 0.0010888991053022614, 0.00204481665392057, 0.0013111327410163204, -0.0015457590510652053, 8.941043612908183E-4, 0.0048989580804333676, -0.0011826748482830069, 0.001522466207862734, -0.0017629091010967462, -4.532324695266493E-4, 0.001354224898038445, -4.4669992070378627E-4, 0.00252033862222889, 0.0021215565517371513, -0.0043846552285104696, 0.0013770297979729053, -0.0033551776661492954, 0.0017920651136101022, -7.670647497274754E-4, 0.003434434693851547, 0.00196464788939608, 0.005141426965432124, -0.008871614241249902, 0.00292536317672306, -0.018608929533038193, -0.02283328128711331, 0.031769963527905384, -0.037030563490911474, 0.01689561491251972, -0.014156160388959251, -0.0015691901669656607, 7.426508679503402E-4, 8.5230025187835E-4, 6.979541807972262E-4, 0.0012448294346978377, 6.612908269642951E-4, 1.6299269591489238E-4, -1.5040265849288902E-4, 5.713875621581104E-4, 0.001627381847953667, 9.425275293556859E-4, 0.0010104618010560174, -0.0013422900055235166, 7.861318325986584E-4, 0.0015034551970783483, 0.0018376378232529082, -0.0030724944723452916, -0.0012965460600625834, -0.00296007117287321, -0.0011883469322517461, -6.982206014918785E-4, 0.0020695536879510117, 9.806819798557625E-4, 0.0018248885315120889, -0.0026922205658910427, 0.0012264837322305626, -0.0019196399337577902, 0.001748935268863814, -0.001944635262467531, 3.319806446436841E-4, 8.606289535679115E-4, 0.004122801209521911, -9.210513425169774E-4, -0.0015487853121110718, -0.0024673792617343517, -0.0011445273083239432, -6.378465360899978E-4, 0.0023885434230045263, 0.0015992118348837122, 4.75878637205771E-4, -0.0025282271753047183, 0.001402946818295407, 2.1391687463588084E-4, 0.0029953649788905844, -3.384566925066949E-4, 9.725992089480977E-5, 6.931551756955241E-4, 0.0049820339303220775, -3.843998364082324E-4, -5.369824074661659E-4, -2.2709615880257892E-4, -1.0929947962990717E-4, 7.095867754516349E-4, 0.0016471010235828238, -6.780302490670639E-4, 1.0211406941395672E-4, -8.24293416492168E-4, 8.995716214474981E-4, 6.376193467806136E-4, 1.219909903795841E-4, -6.715463432022781E-4, 4.4617415452948624E-4, -8.344628979742206E-4, 9.471703834965071E-4, -0.0012628557751814822, -2.102920245981479E-4, -0.003175373831366017, -4.3968727814127566E-5, -4.4223830153753427E-4, 0.00255863238311428, 0.002873041785622483, -9.997552055741538E-4, 7.798835813145914E-5, 0.0057900471362103325, -0.002211618587080998, 0.0010221923833303994, -0.002528720678021205, -0.0014812954014321983, 0.0016573800266357565, 1.8302952427522883E-4, 0.0016641541529183311, 0.003112136045276293, -0.003944487492924882, 0.001862404995786376, -0.003731592521877331, 7.240421846262317E-4, -0.002266383612740257, 0.002508338029140719, 3.544200197486062E-4, 0.005375570330154102, -0.0019563640031286, 2.842640869374068E-4, 7.536538216524518E-5, 8.018817322746781E-4, 0.001217893095297607, -5.137621040248444E-4, -3.0051155248100175E-4, 0.001778766244362489, -7.978628398343408E-5, 0.0011381591419234154, -8.098321429125774E-4, 8.894254174912231E-4, -0.0013534420157725936, 3.390213763359781E-4, -8.766381217859475E-5, 0.002188580557308365, -9.386922929830515E-4, 8.096106827532941E-4, -0.002013453313524628, -0.0018491955729260442, 8.262805820376971E-4, 0.002229973064924653, 0.002139121318723454, 0.001601516186646326, -0.0027832293867512486, 0.0017055713345838488, -0.001099133844265084, 0.00150576045632095, -9.363863197758338E-4, 0.0018038588871016818, 0.0012475321179433257, 0.003631086814775799, -9.023143092655671E-4, -2.1071613529217862E-4, -8.399116839403382E-4, -1.6995769453462132E-4, 5.767219300783582E-4, -3.912761601421915E-5, 4.861750022378455E-4, 0.0012516309460367222, -3.286090402214657E-4, -5.494918627338358E-4, 2.864770511057281E-4, -4.4459407123282154E-4, -6.739277290565462E-4, 3.9543798383121637E-4, 5.674788924019171E-4, 5.207641298229496E-4, -0.0027141215829542997, -1.0530562133940612E-4, -0.0019512176417356476, -0.001321353937693125, 9.747479398383582E-4, 0.0010550318370187493, 0.001545383024728657, 0.0023262240602574535, -0.0028068690407990707, 0.0010507590616353844, -0.0022351528808951726, 6.491318325882816E-4, -0.0014512854732375144, 0.0022544196137935183, 9.79899494493521E-4, 0.002676419315309957, 0.0033868383607684833, 0.009116944696552437, -0.00462441369656861, -0.005755741705558502, -0.0025568070318472426, -4.321696162767376E-5, -0.002762270005811944, -0.001677238257910713, 3.2136168319298617E-4, 9.781682827379292E-4, 0.0015102373272474433, -3.207749961536329E-4, 0.002076220554765556, 0.005219155094770069, -0.00329144156181866, -6.215764609430064E-4, -0.002997113026500205, -0.0017404091449252479, 0.0010752381646047476, -2.7154527735781834E-4, 6.591088055452623E-4, 0.002447592677305246, -0.0026223173545213163, 0.0014140454881305256, -0.003102957862098881, 5.473913126114488E-4, -0.002172842104837704, 0.0022313711980264227, -8.88580037467248E-4, 0.002563819381645136, -9.915497976735663E-4, 7.065575900880843E-4, -6.171568470862495E-4, 5.401417747210409E-4, 3.174454035532919E-4, 7.569664550342639E-4, 0.0019678513943182994, 0.0011484265968219995, 0.001815432545037908, 7.733096768892129E-4, 0.0016768005821313384, 0.0013522895580969092, 0.0011422925424371668, 5.600296569114862E-4, 7.528075102155589E-4, 0.0018835614926896345, -0.0026538522448835687, 2.921805679271584E-5, -0.0034415423259246633, -0.001863304234287367, 9.205949903676046E-4, -0.0020310159536253355, 0.0015399789841121475, 0.002442542884773955, -0.0020170658592301638, 9.761605647042415E-4, -0.0010376340456072332, 2.520394459071989E-4, -7.544438706504281E-4, 0.0017860032076543325, -2.1162609167612906E-4, 0.003983573039493471, -0.002451202135457933, 7.500005658933845E-4, -0.00282866858000258, -0.001978955503128369, -4.409117645219965E-4, 7.376254444700406E-4, 7.755457084210582E-4, 0.0022814419654853552, -0.0028043697997620527, 0.0016661073103447625, -0.0025629831178890753, 9.778147449238019E-5, -2.7586204895180527E-4, 0.0017951584097382965, 6.558920638109441E-4, 0.0034793220423456088, -9.164051075331515E-4, 4.0606742756569593E-4, 4.720937240366653E-4, 6.956663041726209E-5, 0.0012854893861190858, 0.002760783557454837, -4.905147768248367E-4, 7.882326739462731E-4, -5.8360256352142995E-5, -3.558500430227633E-4, 0.00129291128636796, -7.793020373335216E-5, -9.705188685447572E-5, 5.860594446634388E-4, -4.374267249524935E-4, 0.0015267782711006418, -0.0036060607332874715, -0.0018212942126067467, -0.0035821716548993096, -0.0023985059431599447, 5.658401689749484E-4, 0.0011835177002966764, 0.001358406873889502, -1.0922458375556884E-4, 0.00114229886362093, 0.002505077473814292, -0.0023558260829598878, 6.943338982639228E-4, -0.004059083570285307, -0.003537183677395573, -5.321119173130949E-4, -0.001047902410048558, 0.003329888555938661, 0.0021412006821490475, -0.00476506749557337, -1.0097839462955456E-4, -0.002017352266602126, -0.0016699582276526034, -0.0018501221267772019, 0.002806363679589157, 2.532665210409728E-4, 0.0027477681798995313, -0.0014599476577772637, 7.490584245366603E-4, -0.0018399355302749954, -0.0010750724413781487, 1.7481377386878503E-4, 0.0019773060210648207, 6.601674843211619E-4, -1.632474935365789E-4, 0.001002713010773021, -2.446797030807755E-4, 0.0013829401475379862, -7.026515156235199E-4, -7.22964336571695E-4, 9.37108232865378E-4, -0.0010630733445850058, 0.0011304987216271378, -0.0015571487370496758, -5.567712033740448E-4, -0.002705938439689681, -0.0015713472100279643, -1.4887522352338917E-4, -9.477901156616339E-4, 0.0019564087907539556, 0.0029200494686543064, -0.0023379598062549644, 0.0017562451812234044, -0.0019835653516116825, 4.8638500880029045E-4, -0.0015193185020441612, 0.0033229625818430968, 0.002015039299822056, 0.004016207488690448, -1.7978422520293844E-4, 0.0012952828110877226, -0.0010724942983976202, 3.14825483464259E-4, 1.9667128679554525E-4, 0.0018876324534007633, -9.678406760860465E-4, 4.283995369953582E-4, 2.7802923290165028E-5, 8.634361468350788E-5, 7.131550016227593E-4, 9.648728686334907E-5, -1.5198369902413812E-4, 0.0011158252218101674, -2.5577548245029735E-4, 0.0015652905931014633, -0.0028685430134808717, 6.985696821016052E-4, -0.002890335878090955, -2.853248374366166E-4, 0.001939977057991553, 1.6707175394864654E-4, 0.00392384433705841, 0.002766017507162934, -0.0015654473982309799, 0.0017604972643715075, -0.0013548905175567235, 0.0010819916708947233, -4.646397597949585E-4, 0.003037510972196288, 7.346798812031805E-4, 0.0034148308532292626, -0.002042762983264255, -3.0260426468313646E-4, -0.002766703259101299, -0.001162892965390617, 0.0014066779305269748, 0.0014481161764130662, 0.0011587972056448061, 2.279283819718106E-4, 6.755328796896496E-4, 0.005533763847527717, -0.0011431232418033632, 0.0028363728249283645, -0.00273483445533035, -0.0010303497276013584, 0.0019938450319157155, 3.220401653435136E-4, 0.002264027777387638, 0.0032474659695960447, -0.0025194234340082397, 0.0025522107109493056, -0.003305582471911628, 3.023913341601383E-4, -0.0010634397471043706, 0.002739318013766363, -6.417729816888013E-4, 0.0036677554340088904, 3.497934603052073E-4, 0.0013036722186177312, -2.0117246696889184E-4, 3.8175941837913403E-4, -3.658260237937997E-4, 0.0021429665652037844, 0.0013180308800204241, 7.129469966277208E-4, 0.0010885157728513606, 8.493433344169583E-4, 5.706233726825508E-4, -3.9714331370796583E-4, 0.0013123323838678564, -3.546042260748716E-4, -2.2466497205385978E-4, 0.0027952225039056307, -0.0010490786825269186, -1.7004196206142782E-4, -0.002822508569866263, -0.0013760715388381568, -4.1939692207209375E-4, 0.0013908673988578062, 0.0010672796422315006, 0.0022041513413312347, -0.003010146261378751, 0.003037287224282699, -0.0010487547332632798, 2.7017692429917063E-5, -0.0014002535362696907, 0.0022279555293226805, 6.661949509966562E-4, 0.003749322024916535, -0.0025982896386898876, 0.0011249126828266724, -0.002877096070780032, -0.004204821967070062, 7.974872592818542E-4, 0.001079305856126742, 0.0027017580311798793, 0.001909806634217739, -0.003069239238547571, 0.0015623405321721446, -0.00195771532427554, 0.0016489092104287548, -0.0012835859367541684, 0.003142360605779694, 0.0014193018472340503, 0.0031572268619601803, -0.0015394981117135623, 6.583830820951647E-4, -5.004389228558343E-4, 2.876507856829373E-4, -1.395977124590981E-4, 0.0014382581275022672, 3.102350071431991E-4, 3.7463997995169007E-4, -5.370532890003484E-5, -4.168367344854762E-5, 0.0011343354459828828, 6.715538332669338E-4, 8.954603651584868E-4, 0.0011226304828661482, 2.8313709339770443E-4, 0.0019409679777339685, -0.0010319878942326105, -6.577739520280531E-6, -0.0019657892699636217, -0.0019042793992621603, 3.5153505813353926E-4, 3.0322038395703185E-4, 0.0010433593219954943, 0.0031415394950786175, -0.0024291089668444818, 0.0031570663593813346, -0.0016573425513413654, -1.1085740995544612E-4, -0.001892244299526649, 0.002840998901779393, 7.009582500379095E-4, 0.003939651905697103, -0.0015032209840699193, 7.358190493069879E-4, -0.0019073581022477368, -0.0017247220224473813, 0.0011299125508814239, 0.001590203677553293, 0.0031284220410143053, -8.974987836253066E-4, 8.91610426098465E-4, 0.004636037571935151, -0.0015366691146008844, 5.375427984159395E-4, -6.147930107323375E-6, 1.1504429700950631E-4, 7.610659388973193E-5, 6.877829832393385E-4, -1.014118647193695E-4, 4.320197033589469E-4, 2.0726148903821295E-4, -2.952363836438935E-4, -0.0010957799452668187, 9.011874835616169E-4, 8.213072789144613E-4, 0.0011793265194395479, -5.152994112224705E-4, 0.001373769703027104, -0.0019408539135609016, -4.390310987735615E-4, -0.0022310753563632006, -0.0020326538593286457, 7.696188929784139E-4, 0.001388825203062793, 0.0018985398062372873, 0.0020975046473516094, -0.0020115731341142956, 0.0033164597017935828, -0.0010161656367373545, 7.08655342843583E-4, -0.0018160835888874176, 0.002950142408156846, 2.7860563243183815E-4, 0.002027889157433897, -0.0029999661114340154, 1.4736322619954398E-4, -0.002384931133657491, -0.00273765052273975, -9.886034292179956E-7, 0.0013427975831181663, 0.0011086707455408309, 0.0031300640422842634, -0.0032139003089188462, 0.0017974904518578663, -0.003267889757340091, 2.754508980700339E-4, -0.0023131208484391784, 0.0024221469172674027, 5.969520890587231E-4, 0.002043676180320945, -6.879274852957017E-5, 5.565658222651121E-5, 9.74230023290338E-5, 5.523569550932301E-4, 2.855491144688224E-4, 7.059953107468675E-4, -1.9109994095713656E-4, 0.0010123921840098607, -6.735877359234947E-5, -5.397799390500014E-4, 3.7353444503954256E-4, 4.051737939149481E-4, -2.4598376871792804E-4, 5.013841219086623E-4, 2.644657436822802E-4, 0.001257963370546302, -0.0012323945631681084, 4.588165508000461E-4, -0.0031657247638495934, -0.0029523680099904, -2.001099257633617E-4, 2.005109579340126E-4, 0.001161104544276312, 0.002198088361829007, -0.0037039198636482458, 0.001668728472290136, -0.0020271872013493952, 6.790080244498723E-4, -0.00140004212610112, 0.0020598001388657966, 5.60113775746602E-4, 0.0029521943523076024, -7.998917741280302E-4, 0.0010645286723133128, -0.002516316004912424, -0.0010398768276062647, 0.0012358609369712437, 0.0017853170749201931, 0.0019027166298355254, -9.409453592201826E-4, 0.0012682994138181502, 0.002887617720531309, -6.089952549732881E-4, 0.0013067281847675887, -7.407128326247544E-4, 5.480068396382005E-4, 4.677102970727155E-4, -3.059665576936347E-4, -3.5520736011898386E-4, 0.001684221614489955, -9.319273988253047E-4, 5.951384383151406E-4, -7.741908747801781E-4, 6.882838251239874E-4, -6.241170438630963E-4, 0.0016361119456472094, 1.7450439396947837E-4, 0.0016823210237272438, -0.0022699898919904507, 2.61235350472368E-4, -0.00367523956140797, -0.0017421050040962237, -3.390991951057878E-5, -1.7436088017820108E-4, 0.0013188935696842644, 0.0014424377241962767, -0.002725012676217278, 0.0012719114190217558, -0.0012187268392603499, -3.809227127847298E-5, -3.715761483063353E-4, 0.002240363976105408, 6.497024899171084E-4, 0.0035335417365231917, 2.813410495255868E-4, 0.0019193990352872447, 5.730837297900336E-4, 6.217939287165689E-4, 0.002818706453284859, 7.451365094499075E-4, 3.053065104579298E-4, 0.0018767174091494365, 0.0010638558609661795, 0.004269001073191233, -0.001842594709174712, 0.00116913804818496, -0.0033759978372744437, -3.304879288407645E-4, 0.0010972915723890533, 4.2374134086837524E-5, 0.001815757070938663, 0.002386238068186903, -0.00228985751652777, 0.002220004558041502, -0.002074737360145993, 7.933381345309395E-5, -0.002357663393727446, 0.0018580958101608656, 0.0012045317301979258, 0.005532086432333394, -0.0012127387683460843, 0.0019518249233816468, -0.002401598531779009, -8.712862983878836E-4, 9.229424598240724E-4, 0.0029695837506184537, 0.0016949488922092357, 0.002501793506514388, -0.0012199652246011527, 0.001967231786158677, -0.001825437178256745, -6.970739050308041E-5, -0.0028365985262714458, 0.0019394007038442164, -5.382367953106759E-4, 0.0035692013797911197, -0.00162470716402435, 3.292191182071467E-4, -0.0036948594236418195, -0.0010091495896415563, 0.0012164170160328999, 7.110572016166888E-4, 0.0031327261265722615, -0.0017069532670820728, 5.434982949292805E-4, 0.0041443655365129405, -0.0020383841799408143, 5.755630239629734E-4, 3.554873176783409E-5, -7.879206621002775E-4, -5.54203386120222E-4, -6.333291205463248E-5, -1.8820047303714897E-4, 0.001611184205059865, -3.0471187571188673E-4, 2.3336089183600148E-4, 0.0016897226721354452, -2.445472920629653E-4, 4.1738330998491134E-4, 3.521526358262932E-4, 0.0012606276174891953, 0.00308564268254043, -0.0014921788341862078, 7.880470028338214E-4, -0.0021256094919910686, -0.0011216488062306715, 5.634670395442132E-5, 0.001722049045642348, 0.001592557095410869, 0.0023542290222497097, -9.372577647771038E-4, 0.0017284398374480914, -5.669331902117417E-4, -8.660264424897347E-5, -8.38806267312776E-4, 0.0017128085183627535, 5.73879962418357E-4, 0.0036413502680295905, -0.0019526410651862535, 4.022494903516212E-4, -0.0019058528388899711, -0.0019383597340970633, -1.7318453937621128E-4, 5.3682727471844E-4, 0.0023076797968791135, 0.0016110300303559836, -0.002404533981602405, 0.0029915917813358845, -3.785313657907359E-5, 6.821188072759512E-4, -0.002371802100925886, 0.0021096274219502493, 1.7829572429472963E-4, 0.004154744636953619, -0.0013968461979224297, -8.589245059655112E-4, -0.003481790171327819, -0.0019996772364104165, 0.0011219111732680598, 0.001257278418452252, 0.0025521313383531754, -0.001615131338178804, 0.00223565713677201, 0.006247756426780202, 0.001963255223587748, 8.886604034105078E-4, 5.862307008162788E-4, 0.0014295438012048153, 0.004795072135391273, -1.52563241836215E-4, 4.213519887520914E-4, 0.002909277853241547, 0.0010510227678201833, 0.005868424890362559, -0.0015558822917829613, -1.8457845596016147E-4, -0.0028790912229879253, -0.00212761935924558, -5.817969840415469E-4, 0.0017212951869608533, 0.0028745897129176727, -0.002342484242396387, 0.0015889857481514846, 0.003015763361258054, -5.858652529675561E-4, 3.9583504990357604E-4, -0.0019315556353927976, -5.556828939986213E-4, 0.001797868089831367, 0.001503200696718208, 0.0015456630796276383, -9.733679415363889E-4, 8.580560523487402E-4, 0.004558038205276545, 0.034405433229743605, -1.3716188043256584E-4, 2.890306872773931E-4, -3.438221148594921E-4, 3.9890408016320545E-5, 7.421271455449047E-4, 2.0443469112286974E-5, 4.4553416729246507E-4, 3.647415365663266E-4, -5.932240903255029E-4, -9.407748273533643E-5, -1.1772998410292015E-4, -5.1323755392143513E-5, -3.7972583538348675E-4, -1.704674922714347E-4, -6.41111696439651E-5, 1.1421083120641891E-4, 9.190847335065525E-5, 6.307476775673946E-5, 3.8089853613563136E-5, 1.1058855053449744E-5, 2.70828306855855E-4, 2.9664388090392963E-4, 3.5108121246953865E-4, 2.0936278273929386E-4, 3.6791344378784224E-4, 3.0951947407098826E-5, 2.983598972497556E-5, 3.442051601171731E-4, -4.747931154620967E-4, -1.9745997610084895E-4, -2.6786471157293256E-4, -3.694311538740191E-4, -5.314637528130874E-4, 2.552487617234838E-4, -6.08994335786942E-4, -2.221824783480834E-4, 1.1027425845912232E-4, -1.6064239382237934E-4, -1.465053822225191E-4, -4.150934776842984E-5, 0.003303740573576437, -0.0010026775990158597, 0.002546472309693484, 0.0011606856461887068, 3.9641739928669185E-4, 9.969702352584425E-4, -0.0011383690098659904, -0.003186867361093899, 0.00539091071981295, -0.0011386174103695558, 0.005543832319788177, 0.003963514313914767, 0.002003641898445485, -0.0041822265645212135, 0.0033414980593254413, 5.380663899029953E-4, -1.2256357850944912E-4, 7.84126861765632E-4, -4.1451435233956675E-4, -2.11205121252156E-4, 1.0575324271841599E-4, 3.431675191106205E-4, 4.680301490291492E-4, 2.1116833627510056E-4, -1.855642620361746E-4, 2.748969887717801E-4, -4.666204599222115E-4, -3.0409342033748933E-4, -2.976162847458672E-4, -4.516970982868111E-4, -3.029769338863196E-4, -3.8847116149688714E-4, -1.9120220939194786E-4, 5.176195934365698E-4, 5.258977297766489E-5, -5.725322661113385E-4, -5.268014828860296E-5, -1.961516375204046E-4, 1.8963642036474014E-4, 3.9851608072989826E-4, 1.646642771458316E-4, 9.923087237443285E-5, 1.7869754088594336E-4, -9.601603001563319E-6, -2.8146838842317047E-4, -3.5340489943973925E-4, -9.520084869756376E-5, 5.086961468302009E-5, -6.220614322627735E-4, -2.0137231005244597E-4, -3.1606066155797746E-4, -1.3066788300170138E-4, -6.037970622454396E-5, 2.6140177594620884E-4, -3.414417918027882E-4, -4.717434622379191E-5, -2.307244586618415E-4, 3.2361815503395433E-4, -4.841281053596687E-4, -2.934698061720885E-4, 3.6740802380637566E-4, 5.5470536755679355E-5, 2.2644615824816386E-4, 3.238287173865762E-5, -3.195067767569013E-4, 1.4522475031090065E-4, -5.033261412856579E-4, -2.968239074023981E-6, -2.3335934216078287E-4, 2.893823479763245E-4, -1.7093888503729246E-4, 2.085292868203634E-4, 2.533740305345726E-4, 1.5712361129334744E-5, 1.2373465432093876E-4, -5.710211591971431E-5, -3.428892180764348E-4, 4.7575133676900445E-5, -1.743479050830802E-4, 2.1060131398055972E-4, 3.966406634927175E-4, -1.085130193522265E-4, -1.2241782540282876E-4, 3.00737505020735E-5, 2.1212299347794095E-4, 5.795776564591882E-6, 2.1855436273917283E-4, -3.337476234425675E-4, -4.908522472625225E-4, 1.641669463979384E-4, -6.331509341984921E-5, -3.31733049011021E-4, 3.136633420605958E-4, 2.2312819058555492E-4, -1.3656859523735626E-4, -3.0645063019910367E-4, -4.976601152620894E-5, -1.7595765637279902E-4, -2.1950302496667783E-4, -1.4228557246045385E-4, -5.114122783625456E-4, 1.1525936693119909E-4, -1.8327996995628224E-4, -2.460615863099742E-4, 4.359111137044623E-5, 2.327951518452009E-4, -1.6300716876245513E-4, 1.4823471090450236E-4, 9.27011188893433E-5, 3.69780515028155E-4, -4.035370666118414E-4, 1.790832428393977E-4, -7.736756047722215E-5, 1.8706160129425622E-5, -3.4819121203428333E-4, -6.58363892987795E-4, -8.100415240507318E-5, -2.0337401439982647E-5, -8.116461849917527E-5, -2.7908611217392846E-4, -0.0063713328978554015, 0.004711021304980381, -0.006119837205792996, 0.003454628205230705, 1.8625995362000765E-4, -0.0024122576895004143, 0.0051869910910464735, 0.0033954103897275425, -9.130532265661287E-4, -0.0021669353610949105, -9.026975304399727E-4, 0.0023650508963636617, 0.0036431151339099696, 0.003033088671171995, -3.474724479666124E-4, 1.5317132411190825E-4, -2.244734063175846E-4, -5.57704024715783E-4, -2.4129805764809678E-4, 3.589822618187207E-4, -4.160256276186571E-4, -3.471603258706549E-4, -2.212521467111E-4, -9.459570769809636E-5, 2.503925876469582E-4, 2.314301182885825E-5, 2.9191514933313756E-4, 2.5986821588919987E-4, -1.9365954705831303E-4, -3.5495333982226687E-4, 2.508837821465539E-4, -1.0352683080563018E-5, 1.6919277374764763E-5, 1.9928395074401322E-4, -5.057829889959678E-4, 3.450505958066093E-4, -3.819329470532998E-4, 2.0920256992083656E-7, -1.966900828110433E-4, -1.6198401546889342E-4, 1.2167109367023754E-4, -7.45461987142458E-4, -3.22956826806055E-4, -3.401820300878352E-4, 6.274740920482691E-5, 3.3342220430792224E-5, -7.815065595693614E-5, 1.5490806981619714E-5, 2.905531200134393E-4, 2.4451145334648035E-4, -4.4209118802045926E-4, -2.659933886569605E-4, -2.2878058215897667E-4, 2.798176122008561E-4, -7.720148340477774E-5, 9.531756388134519E-5, -2.919002050346285E-4, -8.200001275519977E-5, 3.3231964482454313E-4, -8.897457196505783E-5, 4.30603327098376E-4, 4.725843328128596E-6, -4.903774266314974E-4, 1.547882076548239E-4, -1.5044989493761987E-4, -2.065221034219494E-4, -7.09423551115716E-5, -7.795395247978815E-5, -2.6502859596515613E-4, 1.8538915777760576E-4, -7.090011096063726E-5, -1.7264185822082976E-4, -1.3567565807506019E-5, -7.994520919128186E-5, 1.627946139670685E-4, 1.2138785301151205E-4, -4.046995323885326E-5, -8.343623178873327E-5, -4.690170991094331E-5, -2.4330815382828737E-4, 2.1704739788075628E-4, 7.813234865777993E-4, -6.612358095204896E-4, 4.011857425326431E-6, -4.462947850227997E-4, -2.919440264107357E-4, -3.697960022461668E-4, -1.1003701326497842E-5, -5.145835634182334E-4, -2.467655488226559E-4, -2.4034572429607547E-5, -7.954476022325443E-5, -4.6633428941081833E-4, 1.3199315078679023E-4, 8.413721577084292E-4, 7.075165089839351E-4, -1.7993418472454358E-4, 2.1246550252910666E-4, 6.133434059049573E-4, -1.3281464993693097E-4, 4.8749733746286943E-4, -1.8066089686768832E-4, -1.4883171261777917E-4, 2.433914765803007E-4, -3.700596939259299E-4, -4.391934877167366E-5, -3.906806299402057E-4, -1.2020904463072183E-4, -1.587739749800017E-4, -1.1049793601000674E-4, -5.336639483097434E-6, 1.4128308811887003E-4, 2.969360461080307E-4, 1.0943591682130178E-4, 1.5551280378138836E-4, -3.765044986323853E-4, 3.885463443780697E-4, 3.5440047605291936E-4, 2.1529340544848888E-4, 1.8670576153152866E-5, 2.023769114069539E-4, -6.0396650378150663E-5, 1.3885159317459132E-4, 2.1575552499420597E-4, 5.26390638098251E-5, -8.944107476709226E-5, -6.147743345085129E-4, 3.159843077543194E-4, -2.330748744795168E-4, -3.173441971304374E-4, 2.8655538400438E-4, -2.769401244367279E-4, -1.1869404243988941E-5, 1.4031084985856703E-4, -2.6620807641427604E-4, 1.4219401699020633E-4, -4.2601570627528833E-4, -2.2531330708065787E-4, -1.3940479079045572E-4, 2.470937490672946E-4, -3.650652373939848E-4, -4.030946004183482E-4, -2.758343033850753E-5, 2.934302944837222E-4, -2.6760103827041945E-4, 3.1373455778739383E-4, 2.864816251624139E-5, 6.627813725706021E-4, 8.898904216019584E-5, -2.7262390291348676E-5, -1.8421136927200478E-4, 2.5414305132461336E-4, -1.7801885049196212E-4, -3.7312147532007596E-4, 4.844683669172605E-5, -8.65549760137379E-5, -1.6805086089419978E-4, -6.642339153607926E-4, -4.949093500963148E-4, 3.289690367403321E-4, -3.296342735596862E-4, -1.5051130422975208E-4, 4.750200464185009E-4, -5.891926783701087E-4, 2.4809007806431127E-4, 4.396451636984325E-6, -5.167921010494174E-4, 6.52532842150097E-4, -5.953496154114692E-4, -2.8847589222696495E-4, -6.149789939840108E-4, 1.110413466440947E-4, -2.2272506121697708E-4, 2.3856094230849303E-5, -2.193587082991299E-4, -6.321341227304182E-4, -3.853151972607939E-4, 1.0669593176423334E-4, 2.8265288733501297E-4, -5.806080800268763E-4, -4.1806482766168474E-4, 9.841760923053134E-5, 2.2535392634126622E-4, -2.799443065247629E-6, -1.4537846888378636E-4, -3.3814183859155647E-4, -3.12905844473047E-4, -3.814514830718193E-4, -7.937611787311259E-4, -1.8317174062326832E-4, -1.283025845533821E-5, 3.2394956436370326E-4, -7.688812173890188E-4, -5.518144400134956E-4, -5.096522702219229E-4, -2.6998941654608733E-4, 1.0260641538255549E-4, -3.5236392573833675E-4, -4.736805638185258E-4, -1.017491943168077E-4, -7.714461043800727E-4, -5.61503561857748E-4, 4.925064693941041E-5, -2.3235229558247106E-4, 1.9553120285804648E-4, 2.612980456859736E-4, -2.196799737248303E-4, 4.7237041729163566E-4, -4.6086339654878507E-4, -2.840330688947311E-4, -5.119101866353445E-4, -2.4169526371143121E-4, -3.346805348071892E-4, -2.2871475742361106E-5, -4.718415124187928E-4, 1.0825996710059457E-4, 1.8890771395273752E-4, -5.632118354109703E-4, -3.28729075339892E-4, 1.5200815030156885E-4, -7.645267870658463E-5, -1.1794462788726804E-5, -5.8321706847752814E-5, 3.057580638655015E-4, 8.657559497505989E-5, 3.316564739105402E-4, -6.497925069186085E-4, -1.4709926165875367E-4, -9.323334921701606E-5, 6.970235180106353E-5, -5.149454491085626E-4, -2.2831594139030142E-5, -3.693452013705642E-4, -2.188031505614323E-5, 7.756830189473542E-5, -3.075978818129312E-4, -4.6203145099759854E-4, -5.090711039391088E-4, -3.532386714981031E-4, 4.6841991155308353E-4, -4.7443971066603123E-4, 7.132328701605511E-5, 5.571797502219302E-4, -1.8276016481830323E-4, 2.1704906770388924E-4, 1.7151344504227526E-5, -5.002501928558503E-4, 3.2644815569289024E-4, -4.504160243575013E-4, -1.611551202594572E-4, 5.475205714987016E-5, -1.632799861462053E-5, -1.1237026923742422E-4, -4.460734678617931E-5, 4.6970658908684405E-4, 3.6610135164601057E-4, 4.027776467801837E-4, 1.4722983093889795E-4, -1.461937767652935E-4, -6.492679475830431E-4, -1.429627780267235E-4, 6.091670678880906E-5, 2.2486461949803796E-4, -3.635089518814203E-5, -4.917302575640327E-4, 3.3660275604566284E-4, 2.831406925751596E-4, 6.005396940920351E-4, 1.70492614182843E-4, -3.685907194580476E-5, -3.0382863801151896E-4, 1.5327993979289572E-4, -4.778122875991789E-4, 5.987483793513682E-5, 3.03239326739111E-4, -5.410825930183981E-4, -1.5328103976608865E-4, 1.6831682769710422E-4, -1.5270641771907185E-4, -1.211190089741468E-4, -5.746807883740828E-4, -1.7049135224598546E-4, -3.729141864130225E-4, -1.591403163230705E-4, -3.133893826620307E-4, -1.877890056316072E-5, -1.8705249350829635E-4, -2.2948667540883923E-4, -2.3193594292004435E-4, 1.3079487832266638E-4, 2.3800737177579E-5, 2.7837078099263727E-4, -4.259761268422975E-5, -6.36078317155377E-5, -4.408828447311218E-4, -1.5909811657618074E-4, -8.916310235646996E-4, -6.512382130728894E-4, 3.59950194226972E-4, -1.7540432070243422E-4, -5.231634281539997E-4, -4.547054588378214E-4, -1.6396865184023247E-4, 6.777323285357431E-4, -8.339094513200287E-4, -1.0378126798541923E-4, -8.130024480867558E-5, 7.05341222857452E-5, 1.5913913739164088E-4, -2.9865113476259413E-4, -2.294995855865701E-4, 4.399811649686991E-4, 3.7035274500227223E-4, 3.632244346678999E-4, 3.5442964278575514E-4, 1.6737606617704378E-5, 9.097011548608561E-4, -1.9296891890525133E-4, -2.2126638544620568E-4, -3.2845806663842764E-4, 1.911901613377558E-4, 2.791557608913806E-4, -4.436366510945742E-4, 4.813683390317872E-4, -4.106222956285286E-4, -2.222131021646374E-4, 1.2736258383936987E-4, -1.045107083493999E-4, 8.068811520259745E-5, -5.56782518368293E-4, -6.978095835659013E-4, -1.9802289462528917E-4, 4.550895484939501E-5, -1.9549302867700964E-4, -2.4147916301966326E-4, 2.2560782053989667E-4, 1.7718455392420646E-4, 3.584521813907945E-4, -3.22282613652762E-4, -3.8724772146748337E-4, -1.803341530965723E-4, -4.6562067144179336E-4, 1.874382212651714E-5, 1.6996017570093143E-4, -2.901749615735435E-4, 5.829756668173692E-5, 2.730688278474939E-4, -3.477991512786463E-4, 6.728699783882834E-5, -8.170492656125525E-5, -3.5176576999548933E-4, 3.9712283811001533E-4, -2.3413927421371452E-4, -1.2459124073346253E-4, 1.2884931158354358E-4, -1.0045456969270884E-5, -8.441641027045262E-5, -8.145525413301837E-5, 3.920315728334024E-5, -1.678678496466158E-4, 1.446342286461787E-4, 1.8840158424635837E-4, 1.5299089636489446E-4, -2.283739773574171E-4, 4.3581561312976786E-5, -1.4456004394417497E-4, 2.3469385750877792E-4, 1.0913240150629247E-4, 1.480419015087958E-4, 5.583054262509358E-4, -4.1546405672813634E-4, 9.893233744717804E-5, -3.5910279367622005E-4, -2.3095868177488077E-4, -1.3950257444181E-4, 1.6194827445676144E-4, -4.3869078977050576E-4, -3.4146401354028964E-4, -1.5886200138321923E-4, 6.145405551960916E-5, -3.5516989451920953E-4, -8.809785522380161E-4, -2.271326799191518E-4, 7.042255378717638E-4, -1.0549865634333199E-4, 7.01715422522887E-5, -3.094552975353269E-4, -3.548958387675111E-4, -2.8038489546653613E-5, 2.480494302755925E-4, -8.196224036571773E-4, 5.002952005530817E-4, -2.1892444944669309E-4, -1.9244623578048962E-4, 5.556906419518774E-6, 6.058467271423076E-5, -3.836213715294655E-4, 4.0908452980896676E-4, 3.735857579735763E-5, -4.6736819155612314E-4, 3.782897617974762E-4, 7.437566967831086E-6, -2.095076959872617E-4, 2.9821395574970017E-4, 1.3120514847976493E-4, 1.0408818066829132E-4, 3.7872588196855386E-4, -8.263195765038206E-5, 2.5446636852531764E-4, 3.345421872713094E-4, 1.7877749067100288E-4, 7.472562339750911E-5, -6.037510505068575E-5, 8.853926373054229E-5, -4.209691340486165E-4, 6.326735086215716E-4, -2.930594010752784E-4, -1.1383387009237611E-4, 2.7975374605560335E-4, -1.3206592198309382E-4, 2.6771634725394336E-4, 8.2962780320758E-6, -5.321790647797456E-4, 2.0051978022892901E-4, -6.430321174310865E-4, -5.36850724475549E-4, -2.816414644097691E-4, -5.996178979154826E-5, 8.070883804525686E-5, -3.980277760644129E-4, -2.945770974209515E-4, -9.909972085034215E-6, -1.870902657278399E-4, -4.8194719559380714E-5, 3.046527868124385E-4, -2.8059981920113583E-4, 5.418919979143942E-5, 2.0445606183746592E-4, -5.099622855793701E-4, 7.054569498971886E-5, -2.918228103151313E-4, -4.122012331911647E-4, -1.0784179682054039E-4, -3.2105717634732674E-4, -1.3263884764392957E-4, -2.9520804506782353E-4, 8.742968205223672E-5, 3.9161958534808473E-4, -2.923132435013117E-4, -2.822669544395092E-4, 4.68730629270843E-4, 5.079530081064101E-5, 2.1082665951975126E-4, 2.3618348681257606E-4, -7.31340387665322E-4, 4.8609478158736575E-4, -5.155259614693433E-4, -3.505483915712326E-4, 2.6740925193697505E-4, 3.97344441308068E-4, -5.909104349779872E-4, -6.142311407969148E-5, 3.859172481773332E-5, -2.0692592656018204E-4, -2.1235220795704396E-4, -2.2147013597556293E-4, 3.007164790583938E-4, 2.0688518780026803E-4, -2.301154555019477E-4, 2.702188506609501E-4, 3.5100913268294764E-4, -1.7529456851554095E-4, 6.617371039122511E-5, 2.879189346153139E-4, -5.488258525938903E-4, -6.08868917788127E-5, -5.100337928552598E-4, -3.265370972273878E-4, -3.1731328169096683E-4, 2.651360919626452E-4, -3.234971934312314E-4, -6.735625970644738E-4, -2.0066078419785842E-4, -2.216126361619089E-4, -2.3727003782635472E-4, -4.356871450879095E-4, -3.775956134657476E-4, 4.3648521103170455E-4, -2.7324028430449423E-4, 1.7263212089291112E-4, 1.7873323194728021E-4, 3.012215832721274E-4, 6.197845608234222E-4, -3.1206499094946124E-4, -2.6160680916748557E-4, 2.4070809490695243E-4, 2.672878858786488E-5, -2.3988761527678887E-4, -1.2524918272963228E-4, 1.1068149724145695E-4, -1.1722058314522652E-4, 2.1901292043280006E-4, -9.484545942903624E-5, 4.8477761490595334E-4, 6.621693791593113E-4, -1.7353688327007488E-4, -5.955614865428194E-4, 2.726201310450525E-4, -3.0034091081078927E-4, 1.0800390872447302E-4, 2.3977068061051894E-4, 1.0182830536667367E-4, -4.477836109259331E-5, 4.3655317184862637E-4, -3.22327879036335E-4, -1.638321013971728E-5, -7.96569302712075E-5, 4.4648896224591173E-4, -6.783421409194516E-4, -2.5914177126408467E-4, -5.801707361251107E-4, -1.312641952813605E-4, 3.7416035021141554E-5, 2.703393280572479E-4, -2.3673465692797269E-4, -2.659557644332498E-4, 2.992762588280828E-4, 5.576583999087344E-4, -3.876422416347501E-4, -4.991791987335324E-4, 2.9193474331404116E-4, -1.808774823917462E-4, 7.077455861615031E-5, 9.999046940781708E-5, -7.104818204866558E-4, 4.5979721824867625E-4, -2.627097953503642E-4, -2.2084376446421803E-4, -1.6562425768346185E-4, 4.7112923453249903E-4, -4.5302476713260213E-4, -2.7862160815351816E-4, 2.881175328771586E-4, -6.164031986292213E-5, 9.511547754929622E-5, 2.960052747486053E-4, 3.6049278383502793E-5, 2.220238828111406E-4, -4.460815194379875E-4, -1.7898844602096077E-5, 2.4447647835387483E-4, -1.1947775474973577E-4, -3.834982941306312E-4, -7.057479890904775E-5, 1.640573777956383E-5, 1.4865393778297682E-4, 1.6801590056527175E-4, 6.983050872984493E-5, -4.00781607223379E-4, 4.210903660350743E-4, -9.587505452148546E-5, -6.61512155654351E-5, 2.867761681322224E-4, -2.533852361625461E-4, 8.180829646691834E-7, 2.1385280696043516E-4, -3.277490731969308E-4, -4.3095574555256197E-4, -3.8574236330307454E-4, -7.931364404446374E-4, -4.0996071207192303E-4, 1.8884799988887094E-4, -1.1400641348984896E-4, -3.6378682808931734E-4, 1.741777084119029E-4, 2.113304965616463E-4, -3.121108629889552E-4, 7.19360350389369E-5, 8.63524819682497E-5, 3.041641962311445E-4, 1.2371546906707635E-4, 5.841871480060344E-5, 9.19072957101235E-5, 8.495712002045487E-5, -3.6282498839682785E-4, -2.8303688168341916E-4, 6.600671643625404E-4, 2.0784460523555145E-5, -4.8293354675188254E-4, -4.136488813624848E-4, 2.468882886573743E-4, 7.874594484640637E-4, -6.360473190547598E-4, 1.397859017446476E-4, 2.145865392841717E-4, 4.907970345412394E-4, 2.0564431049847074E-4, -3.7936151098703734E-4, -1.3804726071109613E-4, 3.916143085977539E-5, 5.547889514286777E-4, -1.7418645021531094E-4, 3.633091014835244E-4, 2.167264242700028E-4, 4.737729296516869E-4, -2.982650867635196E-4, 3.005504031884428E-4, -2.644432836764062E-5, 9.731380237378755E-5, 1.7203633891142684E-4, -1.9656651128376502E-4, 4.2330147502540476E-4, -9.533134975058912E-4, -8.498691476643665E-5, -2.483084451606639E-4, -2.119860339636231E-4, 1.919538346076763E-4, -9.378649293567377E-4, -5.812063070452956E-4, -3.0120582897338703E-5, 2.142428444056779E-4, 3.38100000274993E-4, -3.518135857084898E-4, 6.414292655964465E-5, 1.4694197828405065E-4, -2.763198421591581E-4, -3.5944388674014863E-4, -1.7493526939716945E-4, -3.149184675253435E-4, -3.12827819308893E-4, 8.425647607697147E-5, 3.246431337820818E-4, -5.814848396886749E-4, 6.649604424477688E-5, 2.0771932134006727E-4, 3.86215132955665E-5, -5.813295255454123E-5, 2.218702933907038E-4, 2.1348850675028492E-5, 3.094047132315213E-4, -4.03149513118882E-4, 3.578279547812146E-5, -1.4156591991771271E-4, 1.1409911413132602E-4, -2.9339965866614746E-4, -3.3359959447556975E-4, -1.9214772657508172E-4, -1.2757523430648393E-5, -3.61413078489362E-4, -1.7014385513827295E-4, 1.216112768975187E-4, 8.235842040922898E-5, -2.7419681764812346E-5, 1.8091199095709678E-4, 2.804391186152553E-4, 8.883607113897298E-5, 2.3001228520175542E-4, -4.205097458604397E-4, -2.795433304195874E-4, 1.8791196380971804E-5, -6.005614498194186E-4, -1.35347036163245E-4, -1.976592409322056E-4, 2.499299905915956E-5, -9.980443786761052E-4, 7.795835263524804E-5, -9.30556878124127E-5, 3.546819765957608E-5, -2.5093538616351445E-4, -4.5086090796212636E-4, -1.9327741888022472E-4, 2.8853659742135614E-4, -6.643305868661108E-4, -3.6911072880366495E-4, 3.6628917092355255E-4, -5.517679766733767E-4, 2.993282490276818E-4, 3.276165505081515E-5, -1.5646555701969385E-4, 2.7303588320401196E-4, -3.36037863713395E-4, -1.1005942824036615E-4, -1.0288747801967025E-4, 1.0128909021438619E-4, -3.8278304567877086E-4, 2.8026183373601032E-5, 4.714732934196922E-4, -2.2480833052895702E-4, 3.936746616021624E-4, 2.370974560724024E-5, -2.6395532729391654E-4, 1.7499061674269563E-4, -4.6569950984826894E-4, 5.703774544305019E-6, 1.624047411000729E-4, 8.629173037267163E-5, 1.944132316586269E-4, -8.304766640640027E-5, 3.029641345909486E-4, -1.4160801985609462E-4, 1.3651889847084614E-4, 4.165306130216224E-5, -2.0521036144916037E-4, 2.87067437006389E-4, 3.166716164942061E-5, -4.1733962469962994E-4, 3.428289947512822E-4, 4.168540341977017E-5, 6.173548945152959E-5, 2.2809250982700296E-4, -2.5825220708477966E-4, -1.7617299881987972E-4, -2.7574643306585585E-4, -4.0908330732659416E-4, -3.408479494367985E-4, 2.453927499792171E-4, -6.73248552223721E-4, -1.4731274217231127E-4, -5.027876943965175E-4, -2.8788509282401918E-5, -3.9464650976222824E-5, 9.684042909582651E-5, -1.0904863541434658E-5, -2.3268837547212165E-4, 2.8254789241537984E-4, -7.914810991159157E-5, -1.6408147908012827E-4, 3.2390458255569294E-4, -6.125339522750779E-4, -3.5687101291041073E-4, 1.0086580852368084E-4, 2.5197633555424726E-5, -5.655682479593387E-4, -6.95787876951303E-6, -9.759362314492513E-5, 2.1591083996656628E-4, -3.504410741800579E-4, -1.6179433872792845E-4, 3.290718926215999E-4, -4.608985940857433E-4, 1.6190365689929822E-4, 1.6695305397045653E-4, -2.871604444441714E-4, 3.6159927751784024E-4, -3.824509094816881E-4, -2.100368885848264E-4, -1.7425823690764576E-4, 6.261449388610498E-5, -4.177064077527653E-4, -2.3919189438624587E-4, -6.493225537536394E-5, 4.211167984610148E-4, -1.7994077477399514E-4, -1.0834997599461091E-4, 3.3890325675565807E-4, -1.6644022102455983E-4, -1.2935745354499125E-4, 8.792455847471219E-5, -5.376799400943523E-5, -3.966578939624838E-4, 1.841814990169543E-4, 3.0834946220340345E-5, -2.5853400685759136E-4, -1.4943255362568985E-4, -1.0792820112725914E-6, -6.152151542918693E-4, -5.421329258544525E-4, -1.2470424097004693E-4, -2.68516015199392E-4, -6.074876764151954E-4, -1.2800978206866027E-5, -1.6912902777703812E-4, -6.079524073067464E-4, -7.527842549161045E-4, -2.7937088802712886E-4, -4.009304362563358E-4, 9.133565897601276E-6, -2.5515875983513303E-4, 1.1630760509123528E-4, -1.801666205061953E-4, 2.6066468759500434E-4, 1.0250435606139912E-4, -4.1087820276674035E-4, 2.464439286144035E-4, -8.362895588681142E-5, -5.842703582443449E-5, -8.66452801251482E-5, 2.2083513246494023E-4, -3.4611574313840886E-4, 3.1214619696521586E-4, -3.190022680877461E-4, 3.802019779466701E-4, 1.5779397990070608E-4, 8.200463238443767E-4, 2.7792526570507236E-4, 1.0609645288875262E-4, -1.1358831807070793E-4, -2.2794849337040256E-4, 1.4328950742056488E-4, -4.988806864075257E-4, -4.218022425625406E-5, -8.100818237751471E-5, -4.931074516576442E-4, -2.5025330976494913E-4, -3.9041511718724414E-4, 1.7159765537881014E-4, 1.4518033410779174E-5, -9.563522820549676E-4, -6.024936317904931E-4, -4.2715059453258446E-4, -2.1207452878486632E-4, 8.023877685112726E-5, -3.1169488195170996E-4, 7.724997669623121E-5, 2.1521420177700182E-5, 4.403921793485132E-5, -3.7087434353026455E-4, -1.588276230723571E-4, -1.348134514872327E-5, -5.560446992207566E-4, 1.87780154038552E-4, 2.16955464354409E-4, -4.871201753305866E-4, 5.846589138632321E-5, -1.3828088333734387E-4, -4.99995420329355E-5, -2.595698523582637E-4, 1.7850011698170126E-4, -2.2294943031068081E-4, -4.972514845827664E-4, 3.893320643719705E-4, -2.320702826971898E-4, 2.4177289164043904E-4, 1.9447443950679495E-4, 2.1800906349672146E-4, -5.071515279461913E-4, 1.678988973188287E-4, -5.636165444810042E-5, 6.810360261965992E-5, -1.5706257104188437E-4, -2.9470135709520884E-4, -1.627089770884797E-4, -2.93799128249942E-5, -4.881901306048094E-5, 3.7344811210465034E-4, 2.0068680080005284E-4, -8.050039960711864E-4, 3.284083665089138E-4, -9.433225492055529E-5, -3.2441944025633947E-4, -2.3335285511907097E-4, -6.099438015506897E-5, -2.058243100054182E-4, -2.7451443005743477E-5, -2.9085764469489458E-5, -2.626508781794167E-4, -3.6858700933787094E-4, -3.6668569753672825E-4, -4.1348143983262577E-4, -6.960014575669552E-5, -6.629936282609382E-4, -1.6768023966341972E-4, 9.796242514737764E-5, 3.003106995144964E-4, -2.9974611123542174E-4, -1.0543600908063087E-4, 3.384832374920388E-4, 2.813674802595945E-4, 3.950995613934028E-4, 4.8953525478106385E-5, -3.1784726802593164E-4, 2.7807086999820786E-5, -3.946317886722968E-4, -3.7999282946722874E-4, -3.087606284793838E-5, -3.029475335347737E-4, -5.933763789552728E-4, -8.063905703578232E-4, -4.150198500642002E-4, -2.0615552655703384E-4, -9.662303689124668E-4, -1.9036087262698235E-4, -4.242995393365902E-4, 2.59327467064377E-4, -1.4880522051943334E-5, -7.059582828580062E-4, -5.870137383065298E-4, -1.9370542724855431E-4, 1.1441430325563693E-4, -4.8389760886908975E-4, 4.2287307883385404E-4, 3.8786812654933615E-4, -7.395414188178517E-5, 1.508772813733324E-4, 6.647533918787259E-5, -2.251739905607676E-4, -1.5812764200355618E-4, -1.592055712791897E-4, -5.273229634580761E-4, 2.0201165644390888E-4, -7.027399570267071E-4, -3.0122649576930796E-4, -1.990971282332783E-4, -1.1263235412777156E-4, 1.7255248326248897E-5, -4.3919810964376024E-4, -8.828390559196268E-4, -2.2880949807624806E-4, -1.2741522974926939E-4, 2.707653116437064E-4, -2.943995136870447E-4, -2.03531027959149E-4, -3.6414559427767624E-4, -4.747114640998168E-5, 1.0262558525091334E-5, -4.000952969476397E-4, 2.1261526019961965E-4, -9.841356046890056E-5, 2.0282187694461204E-4, -3.0254930805620147E-5, -4.5102522636956575E-4, -2.0537074302983123E-4, -2.60174140737173E-4, -4.6978456470125115E-4, 1.5244143959157854E-4, 2.141832515853334E-4, -1.5698991027479387E-4, 4.898915709509782E-5, -4.6202506404285297E-4, -5.246757633081357E-4, -1.8812587737294698E-4, -2.2918153001159875E-4, -4.2351474164816736E-4, -1.5971911746170291E-4, -5.74882660696835E-5, -2.7207277847098607E-4, 5.7452505479144806E-5, 1.7067507577106327E-5, 2.4127443585308386E-4, -2.500095649938159E-4, -2.378632985109655E-4, -9.811244937344213E-5, 1.5060924890015153E-4, 2.160240075124825E-4, 7.582858564858974E-6, 1.6228998582248281E-4, -2.260779252292059E-4, 2.3388392456478213E-4, 3.0586451887472735E-6, -5.37480587551735E-4, -3.6443659813610205E-4, -3.9456863778442297E-5, -6.200797426674275E-4, -3.965055373463186E-4, 1.460927488256404E-5, 2.1838013527549893E-5, -6.775584943592423E-5, 3.4247894515460075E-4, 9.268702556086803E-5, -1.5152616983124252E-4, -3.6127700930665887E-4, -2.9229182961924794E-4, -1.6551848221511942E-4, -6.938879748083834E-4, 3.3757982703783184E-4, -3.2288336544322575E-5, -2.5826048912236746E-4, 3.2644723914340506E-4, -1.3499234472485737E-4, -1.529617362139794E-4, -3.7115728349594216E-4, -3.2269382812276034E-4, -3.089509160679304E-4, -2.295177222507632E-4, -2.0737191855200454E-4, -1.8996577248683353E-4, 3.47911637076851E-5, 1.5496690541587984E-4, 1.221593226215928E-5, 2.4175565741307034E-4, 1.021836638546769E-4, -7.383462171559703E-5, -1.4488999437297842E-5, 6.366299370182344E-6, -6.95148816683584E-5, -8.259895706285288E-5, -9.460281129157667E-6, -1.3886327987816147E-4, 2.842340748348218E-4, 2.3090589493762607E-4, -2.533045334887815E-4, 1.0252925569624144E-4, -2.422306194406435E-4, -4.0166337526304104E-4, 2.0820904153770383E-4, 9.788338780100008E-5, -2.1501891168798312E-4, 2.1477060207602576E-4, -3.6847858945667677E-4, 3.5272239728669337E-4, -5.587630800331935E-4, -2.691333348927151E-4, -3.460889906990969E-4, 1.8491508754401567E-4, -1.1406878944522425E-4, -1.0044494848438553E-5, -1.2161244860477442E-4, -3.231556415703304E-6, -2.845863461029635E-5, -8.371977462063316E-5, -1.3085519896836095E-4, -1.7488319123556186E-4, 3.50153268483546E-4, -5.169665678318456E-4, -2.6418273596093695E-4, 3.7782994554317196E-5, -3.506651129599079E-4, -7.36842971860868E-5, -2.9229965292777087E-4, 2.9906665301626494E-4, -3.616859308894453E-4, -3.765987622555522E-4, 3.649038822313564E-4, 5.1794522535228844E-5, -1.2072743948635144E-4, -7.966066905666598E-5, -3.10946245840628E-4, -4.34297404564854E-4, 1.012394090757129E-4, 1.766724255575079E-4, -4.956995803204563E-4, 3.123595098153924E-4, -7.048661991940436E-6, -4.5022249411653405E-4, -1.829948854651282E-4, 2.4389792698092588E-4, -4.937370555575473E-4, -6.22124777308832E-5, -1.1549646500516098E-4, 2.0097753410076795E-4, -4.452039684932713E-4, 6.48580041765282E-5, 2.786538068326269E-4, -4.7529694011270323E-4, 1.1490249383107971E-4, 3.638394091072843E-4, 9.682598617456858E-5, 3.199907562065796E-4, 1.93085245181046E-4, -4.822770188218486E-4, -2.010711019633788E-4, -1.2001959814848296E-5, 2.5904419547870608E-5, -5.637724989769519E-4, -2.4105818964307029E-4, 4.1941328083129256E-4, -6.628936635914777E-4, -3.7792697577889063E-4, 4.3366430986091736E-5, -2.676390508658951E-4, -3.0194606865607105E-4, 2.8755037223676487E-4, -5.321259276413791E-4, -8.293994998893347E-5, -2.529312624604725E-4, -2.0616546828615976E-4, 2.5180565854384827E-4, -6.174955313825456E-4, 2.4579961624184195E-4, 6.599723026216669E-5, -4.1307957658586695E-4, 4.1882532640944234E-4, -5.964018324105383E-5, -3.397480370718617E-4, 2.4870371467353755E-5, 4.850418223515513E-5, -2.4297704726128073E-4, -1.377999773022004E-4, -2.072532965545602E-4, -7.366245038600662E-5, -2.5563896389235448E-5, -2.8776715982805833E-5, -4.172079211210148E-4, -7.015679292596006E-5, -1.231356428566665E-4, 6.46474470528203E-5, 1.4928794503051048E-4, -2.1555333568913545E-4, -1.5561692868365902E-4, 3.880163932716347E-5, -1.380292919812054E-4, -3.4831690274276085E-5, -1.0704278584043345E-4, -2.739805978166142E-4, -1.7456806866860493E-4, 2.488946793582458E-4, -4.3880190719980925E-4, -1.3717576149178154E-4, 1.5446844862408281E-4, -6.635528512681835E-5, 3.0451546171213206E-6, 2.609886185204971E-5, -4.28353326919383E-4, -1.9218334926554364E-5, -6.291004720517066E-5, -3.9873162045637415E-4, 4.0122577184215536E-5, 2.661398674211653E-4, -4.304896967067312E-4, -2.3908634741562265E-4, 1.353167819607246E-4, -2.0428520977570078E-4, -2.203973148894646E-4, -7.032920475580882E-5, -1.2160433442740204E-4, 1.98085301271134E-4, 3.543385833515951E-6, -5.142490351347996E-5, -2.6456586031618084E-4, -3.302997967962303E-5, -4.1845159285176494E-4, -4.145238452769865E-4, 2.9157906998975806E-4, -1.849324674597537E-4, -3.9181221629544776E-4, 2.1825082006921167E-4, 4.812054478250667E-5, 5.114158668562128E-4, -5.067201922890327E-4, -5.796858498447948E-4, -5.2702687874460366E-5, -1.9415424494513328E-4, 4.242596485354594E-4, -1.1404822522101697E-4, -2.0001327721445085E-4, 2.5852127250644953E-4, -1.6817935466434896E-4, -3.2328544161099623E-4, -1.0414776976570299E-4, -2.8439967356600317E-5, -6.329299173589507E-5, 3.534320218402415E-4, -7.675161631851396E-5, 1.535427157782091E-4, 2.6118869198378004E-4, 4.851224197773077E-4, -3.9971740985032106E-4, -2.474042629781405E-4, -1.7662443141834504E-4, -1.0609144636559103E-4, 2.135167711036199E-4, 2.243024344525266E-4, -9.594423127705842E-5, -5.921783789759752E-5, -3.669166769040531E-4, -4.48587535131183E-4, 2.2951164802953298E-4, 1.6203664325929647E-4, 1.303753997540253E-5, -4.942935552148322E-4, -2.941329073182595E-4, -3.608414201249613E-4, -2.0192609032364218E-4, 1.2726738369489782E-4, -2.565378523724245E-4, -3.8275534624665855E-4, -4.2364810642283256E-4, 1.1210676014826175E-4, -8.420976300353075E-4, -2.4032034348887008E-4, -5.568580983864846E-4, 1.5970094103383704E-4, 2.7754803829796622E-5, -3.535387289647202E-4, -2.962144417408939E-4, 6.164595265551994E-5, 9.209526719274462E-5, -2.669262722314889E-6, 2.3983752713728257E-4, 1.006541920961376E-4, 1.7261035252820022E-4, -3.4587982028908156E-4, -5.683685971845287E-4, 3.0970462921537577E-4, 8.819184576170204E-5, 2.703496431715011E-4, -6.50275958635491E-4, -3.566568474650482E-5, -8.720190172887433E-4, -2.4308986875096406E-4, -1.5931703132208833E-4, -6.455704997052655E-4, 1.9013928164430063E-4, -4.8708991487524523E-4, -0.0010828158197630795, -6.451804485502631E-4, -2.579862220303744E-4, 1.3798812164910585E-4, -4.844950271238546E-4, -3.9186889051484174E-4, -2.877444086454802E-4, 1.7364804660541624E-4, 2.5217097811251407E-5, -5.262203453623911E-4, -6.048246945080685E-4, -6.995002225284358E-4, 7.198477223196303E-5, 5.095548558590543E-5, -2.3102040578802058E-4, 2.843684304290752E-5, 2.878714417028558E-4, -2.2929263986961565E-4, 1.8401238761045274E-4, 1.6578502131815384E-4, -6.763897560235009E-4, 2.2916107317056387E-4, -5.043690585340328E-4, 3.933933235051093E-5, -3.9946460365922784E-4, 1.140308877507363E-4, -4.309724109836891E-4, 9.358489207787009E-5, -2.880190822382125E-4, -6.132497102949897E-4, -2.109328380006161E-4, 5.265608120547565E-4, 1.5839655217565534E-4, -3.070037519569068E-4, -3.773544694982013E-4, 1.6689701633410513E-5, 8.8307441404411E-5, 4.151752764681711E-5, 4.663059603727595E-4, 6.182722154339623E-5, -5.755473556660311E-4, -1.233495911076101E-4, -9.951697040140913E-4, -3.261703608089258E-4, -2.3413708792205281E-4, 3.524182894038385E-5, -5.932430168420987E-4, -6.599873320759777E-5, 7.728966641977063E-5, 2.0061228475471218E-4, -2.5195060856212563E-4, -4.663212134887809E-4, 1.566882299946716E-4, -1.661560053813851E-4, -6.58430991961944E-4, -2.1217588299033124E-4, 1.3191473361498474E-4, -1.2641453596210216E-4, 1.833555187743185E-5, -5.11454060871894E-5, -3.6305535581295423E-4, 5.89490127264426E-5, -1.8096170942128196E-4, -4.7286191002110295E-5, 4.08237931864995E-5, 4.9107385970736006E-5, -1.0150575611164405E-4, -3.135566431223852E-5, 6.641906028635627E-5, -5.385275020583004E-4, 4.13387450615676E-4, 6.502160609488594E-5, -3.4864067606649066E-4, 6.658590681094514E-4, 4.2851007793161434E-5, -9.42746757024174E-5, 3.2836357722284084E-4, -2.823523003673916E-5, -1.418568640423985E-4, 1.266510108616886E-4, 2.6394464491427216E-4, -1.5681701211965544E-4, 9.541000870064457E-5, 4.7162862273082804E-4, -7.403725510826754E-4, 2.3847503250119665E-4, -1.8764880475121727E-4, -1.4703969771530438E-4, 3.565806182352599E-4, -1.343452366270218E-4, 9.512470584056921E-5, 1.0553120083998617E-4, -2.509249433720325E-4, 4.794828702428043E-4, -3.5811744327612395E-4, -3.6771787954188186E-4, -1.4049349752562367E-4, 2.051316533675958E-4, -2.7514579351568423E-4, -9.14132641223511E-4, -2.505423200609806E-4, 6.821282223258237E-5, -2.0927199977676628E-4, 7.946793715270185E-5, -1.9797734897411795E-4, -9.956205045066684E-6, -1.319344545964044E-4, -9.175543108918215E-5, -4.110281723076617E-4, -3.216909354865561E-5, -3.7691529114630934E-4, -2.3407535764273524E-4, -1.481737306869493E-4, 2.804695895869513E-4, -4.752459232660819E-5, 7.299671873746763E-5, 1.2324732436602855E-4, 7.740803318887369E-5, -1.6716543385440936E-4, -4.785686335395179E-5, -1.9958835268807287E-4, -2.4857504953878693E-4, 2.0850719237754372E-4, 1.9898714140829748E-4, -5.650988282845494E-4, 1.933223954763336E-4, -3.800631472178019E-4, -2.698950789574297E-4, -2.0614004043839726E-4, 3.2819824042540596E-5, -2.808443680316079E-4, -2.1822519685138012E-4, 8.841493208212273E-5, -3.687734108815681E-4, 6.389318663801579E-5, 1.2835465224857813E-4, 5.618783762337435E-4, -6.368266610726003E-4, 2.3451250425909816E-4, -3.105341375875492E-5, 8.648671315092491E-5, 3.495074613401768E-4, 5.640535457390186E-4, -4.2524755299955556E-5, -8.546251788507857E-4, -8.855371726637792E-5, -6.695869981684179E-4, -3.1092033628837726E-4, 1.2189865215031856E-4, -5.396174869863829E-5, -8.505853883750613E-4, -3.293647951137027E-4, -8.60125137706076E-5, -1.1208608089273778E-5, -4.924747645872116E-4, 8.692036545306233E-5, -9.07244340299076E-4, 3.023880585402499E-4, -6.430934065925086E-4, 1.6560335518959828E-4, -9.747603529286681E-5, -2.323006703450778E-4, 2.54457583285496E-4, 9.640825559371351E-5, 6.42352935181906E-5, 3.3075690800453715E-4, -4.284089672335361E-4, -4.0065464775174256E-4, -3.009187535324683E-6, -4.408138822269877E-5, -4.6125833822007736E-4, 3.812587072177741E-5, -4.3311978638940806E-4, 1.955622706114254E-4, 2.5251996836225495E-4, -1.291634018227353E-4, 9.020759646381695E-5, 3.490968934005761E-4, -4.679711130601197E-4, 4.0082373622787605E-4, -4.431875749279896E-5, -2.383715026014983E-5, 4.5631068447352725E-4, 3.820702105636584E-4, -5.737843769787775E-4, -1.3394204098421771E-5, 7.520349670375617E-5, 1.7314657429323536E-4, -5.005616546397535E-4, -2.850159902485813E-4, -3.241530468858994E-4, -6.18716006674546E-4, 9.07975531058389E-5, 1.1333508414832936E-4, -5.159756939019263E-4, -5.357639577556609E-4, -3.879712665902825E-4, 1.872521883519774E-4, -7.152449816986474E-5, -2.0514336593423064E-4, 1.3143112651529433E-4, -1.255671918923186E-4, 7.681754239964573E-5, -2.7608818760662714E-5, -2.9227209874983773E-4, 2.948991036557501E-4, -6.816855303309438E-5, -3.8374201722747104E-4, 2.1475119971840588E-4, -5.004014621988326E-5, -3.3318537827494577E-4, -2.838847618775711E-4, -1.1504677691139649E-4, 2.1841324142750688E-4, 1.3755570261554548E-4, -2.491647741372936E-6, -2.130727957575997E-4, 8.427629802245816E-6, 2.187514841677364E-5, -2.8886436434684046E-4, 5.190029943983039E-4, -6.88512847555486E-5, -1.2125582776339571E-4, 1.6434630150635406E-5, 2.9258382655716105E-4, -7.742081518422283E-5, 4.347899528000446E-4, 9.13282171659301E-6, -7.036802525720291E-4, 1.5318010033887758E-4, -1.3734319037420843E-4, -1.2852369196679494E-4, 2.24941165653622E-4, -2.2959702205393916E-4, -2.9133946223625383E-4, 4.0772857339293914E-5, -4.980894100037965E-4, -2.751857070333911E-4, -1.9573180800993694E-5, -3.221469962036525E-4, -3.582667120046245E-4, 1.1636448898863722E-4, -3.357222029113183E-4, -2.58784087017064E-4, 2.984120323561748E-4, 3.149301930262071E-4, 3.0283150198661676E-5, 4.498086098652658E-4, 2.73188746704958E-4, 3.676719645517067E-4, 2.455774597456106E-4, -2.6071165589983033E-4, 2.5506926893845107E-6, -2.800438557590951E-5, -1.5330869161136811E-4, -3.636213798916704E-4, -1.1857711060904983E-4, 4.3344418484145145E-4, -3.867117769257673E-4, 3.794224586783166E-4, -2.978071918209236E-4, 3.761055026140431E-4, -3.916168062407658E-4, -9.004020387151786E-5, -1.4363169858040974E-5, 3.232364324854302E-4, -2.288155189120528E-4, -6.078008263082928E-4, -9.190760342621786E-5, 1.7099181214812572E-4, -2.326802084693481E-4, -1.0707755712567157E-4, 1.3887189057299354E-4, 3.1297266568750036E-4, -2.339124559981281E-4, -2.7509240462091024E-4, 8.514244800448363E-5, -1.6842256311280074E-5, 4.2345722695514166E-5, -1.2375298495623645E-4, -6.5564667168796E-4, 3.134632447340306E-4, -5.967629089915943E-4, -1.7995202357498733E-4, -4.537294449517429E-4, 2.0636322089579722E-4, 1.4956353158019969E-5, -0.0012107446828480317, -6.180075932601385E-4, -3.291637710324965E-4, -1.0119592016160132E-4, -1.403466533974378E-4, -2.4589240827243316E-4, 5.207297801855908E-4, -1.2323971473484458E-4, 9.726505168157436E-5, 7.590541041962621E-5, -5.502777749257491E-4, -7.586916611591341E-4, -4.078338880619515E-4, -1.8923185422257483E-4, 2.4351076514264998E-4, -5.217970225204472E-4, -1.1788160730180213E-5, 8.050663547990214E-5, -2.3051345504642644E-5, 6.233054350648565E-4, -2.1004766219087883E-4, -3.08984082060936E-4, 1.5554597509691932E-4, -8.55780888136936E-5, -2.399150135899937E-4, 5.343075261111246E-4, 2.3329587349015012E-5, -8.409803501805416E-5, -3.0852183559288386E-4, -2.7422090307995813E-4, 2.938516716856973E-4, 3.0024560678297433E-5, 2.526633151607728E-4, 4.542563301679458E-4, -7.034844263493155E-5, 2.3389111209763783E-4, 1.292236304902358E-4, 1.3908936798561378E-4, 1.7295148777632646E-4, -8.746540949150299E-5, 6.004810470765331E-4, -7.130931825740228E-4, -3.310493618511355E-4, -3.653469446995546E-4, -4.1120547551317044E-4, -1.6781018957713004E-7, 6.632718988592391E-5, -5.436339803569138E-4, -3.6829022934962383E-4, -7.099530162818371E-5, 1.1381362424730342E-4, -3.616807064365181E-4, 9.689755971175424E-6, -0.0010044624352399383, 3.060209816709252E-4, -4.2060323306629594E-4, 1.8003746975153354E-4, -1.9007863267865208E-5, -4.486165876696538E-4, -5.497426394520348E-5, 3.7886616744008215E-4, -3.290896510100268E-4, 1.8652010857642977E-4, -5.710316675692053E-4, -6.795023878260415E-5, -1.259818546694305E-4, 7.0966944899715815E-6, -4.3822276914104945E-4, -1.0757561812654241E-4, 3.2219695648492335E-4, -2.4753077130768497E-4, 3.158319513567037E-4, -1.3284847350817096E-4, -1.62858182911615E-4, 2.1191614167353622E-4, -2.3899016656413674E-4, 1.8440453343866433E-4, 3.7037121041361425E-4, -4.067507886570493E-5, 1.5744631647850095E-4, 1.4197818081029934E-5, -1.0988415753309676E-4, -2.982177205709159E-4, -6.390212493866996E-5, -8.867864004722759E-5, -3.9506717302587664E-4, 2.7058968320907244E-4, -4.779142100738184E-4, -1.2453854223364578E-4, 3.092001092739864E-4, -1.8088612534742428E-4, -2.0115773912403746E-4, 5.704241924464287E-5, -3.8654070663757484E-4, -3.619071324217674E-4, -8.083206263626123E-4, -3.857687503909078E-4, -6.423618118738591E-5, 1.821830379944749E-4, -6.859670922522566E-4, -5.328553825389826E-4, 1.7487204555404292E-4, 1.524740659692696E-4, -4.644559681193187E-4, 3.1469318208885176E-4, -2.92497845840909E-4, 7.421065945207542E-5, 6.8452790997875585E-6, -1.0204015871484256E-4, -2.3237893109359298E-4, -9.267189205754508E-5, -2.293132507980313E-4, -4.869388599181121E-4, 8.205656735866209E-5, -5.895000314384206E-5, -2.6617347767572783E-4, -3.222524924885056E-4, -0.0010219330798016426, 1.1549713058640865E-4, -8.143682816102594E-4, -8.998118411883927E-6, -1.1388477582424195E-4, -6.188130568056771E-4, 1.8956917429764457E-4, -8.062346947371635E-6, -4.040701659882174E-4, 1.1298932066932631E-4, -6.275325608985122E-4, -7.176845997417198E-4, 5.430034393182732E-5, 2.1557150253117473E-4, -3.900045276372694E-4, -2.861431426070132E-4, -2.110841145783768E-5, -6.983881697576994E-4, -1.8972672714238025E-4, -3.4425089962243746E-5, 1.9724033741993195E-4, -2.5053149127362683E-4, -3.135809514673455E-4, 1.609981861230138E-4, 2.1045147401458882E-4, 2.2559473691831323E-5, 8.909276434127043E-5, -1.2893571552390765E-4, -6.241198135998128E-4, 8.67626227861364E-5, -0.001067795085808762, -8.243355985167952E-4, -3.217658818605057E-4, -2.4372610466410615E-4, -3.153495215241958E-4, -2.474161546458337E-4, -1.6599453715277586E-4, -3.6702278222344796E-5, -6.207506745087054E-4, -6.552393526192161E-4, -8.433802990002083E-4, -5.580996198848024E-4, -7.068986889198896E-4, -1.3080423726695802E-4, -2.9269107879398163E-4, -5.048704176296049E-4, 2.62958587122568E-4, 1.5134421677075183E-5, -4.81827357618214E-4, 2.4391323961544566E-4, -5.723970975662685E-4, -4.343989824829675E-5, 4.4751477512411494E-5, -9.170262145126313E-5, -4.438541348497204E-5, -1.8133023298085612E-5, -1.7444915119694702E-4, -6.352804849927829E-4, 1.708158935952364E-4, -1.3840910332792568E-4, -3.9257226555389495E-4, -7.726377932311507E-4, -1.3264531288576014E-4, -2.3574400448363295E-4, 1.6541014364840941E-4, -7.157383431410035E-5, -2.3291821542360065E-4, -1.249707034498748E-4, -5.233887283900883E-5, 6.0394180154088335E-5, -1.1923947891362733E-4, -2.843267929092358E-4, -0.0010288526070554745, 5.4564016283229E-5, -3.555072921921181E-5, -4.2255204760226423E-4, 9.485522070768203E-5, -6.004146498469814E-6, -9.433909908861072E-5, 5.98485359307242E-5, -5.892720465119692E-4, 2.1178907763321327E-4, -5.009438729677426E-4, -4.1637947189424604E-4, -2.1571585302671763E-4, 7.665608549431823E-5, -4.440836137171374E-4, -4.0392317362028E-4, 5.405055738163839E-5, 3.56224110261763E-5, -1.8982787133191986E-4, -3.236176898950393E-4, -2.1922381408288156E-6, 2.971618390285415E-4, 4.5864838746501784E-5, 2.4827952389451914E-5, -1.7032710719387266E-4, 8.995653438061672E-5, -4.950433401920478E-4, -4.350758015167128E-4, -2.5253350448769385E-4, -4.497335034609898E-4, -3.5908131412679416E-4, -4.1576636354138485E-4, -7.149815525161866E-4, 1.929111002918352E-5, -7.729482988595424E-4, 2.7076448751445686E-4, 2.432162343668241E-4, -6.326214032530644E-4, 9.75231667530792E-5, -8.213920224920616E-5, -4.398948164981867E-4, 8.547175712714423E-5, -5.705313538972881E-4, -4.7479591084329605E-4, 8.224541565265303E-5, 3.21398191388037E-5, -3.880319365514197E-4, -2.4843839170991957E-4, 1.2901776288292452E-4, 5.07094282190888E-4, 7.638838297395489E-5, -2.937062346559931E-4, 4.203298148170191E-5, 6.286940744409703E-4, 3.972984817197521E-5, 3.0120043791917727E-5, 2.539947386277549E-4, -4.030132232816347E-4, 2.8342511506372842E-5, 2.8988587779109776E-4, -7.748596800531285E-4, -9.475795512808578E-5, -1.4253507890586925E-4, -3.909728141356778E-5, -6.248579696051294E-4, -1.9457992171488238E-4, -5.446717309107751E-4, -2.8380326968931104E-4, -3.1004992159888804E-4, -9.637957045901526E-5, -3.7846320342065526E-4, 3.0938106886655716E-4, -7.187851013426227E-4, -2.6443488897461703E-5, -3.2455814301965874E-4, -5.1951507127836555E-5, -1.23172047286905E-4, -6.465273888978419E-4, 5.168461603039245E-4, -2.846120294863546E-5, -8.447665677178383E-4, 3.4729843795856753E-4, -5.676654724080859E-4, -3.7634386689636205E-4, 3.353283057338245E-4, -4.427472402910725E-5, -2.9804847815341994E-4, 2.5718174496933106E-4, -4.7942796469367443E-4, 1.985017487066258E-4, -2.377315924488088E-4, 1.962057373391341E-4, -1.0778078140230931E-4, -3.4996253781429103E-4, 3.8297277708933773E-4, 4.760305449883765E-4, 1.6005752854934705E-4, -2.0157363224721392E-4, 3.4551085535890625E-4, 2.3245470040283859E-4, -5.654352485192794E-4, -8.964377683077975E-5, -5.295849209513838E-5, -5.445095118493839E-4, -2.362386559695462E-4, -1.740969967537318E-4, -7.530462393484E-4, -4.578898984611856E-4, -7.208410229217672E-5, -1.5233881915526789E-4, -2.7151334852407963E-4, -5.516873845389806E-4, 4.695854179584375E-5, -4.911433939528042E-4, -5.934653701435658E-4, 3.219475892115317E-4, 3.5923703058468376E-4, -4.87858385264478E-4, 3.7058461173345396E-4, 9.850861250697632E-5, -6.043725165056696E-4, 1.1319668742250177E-4, 9.021732121147855E-5, -3.2758437082986536E-4, -7.427047515234171E-5, -2.039133975997008E-5, 1.7839496963466478E-4, -1.8754395745822238E-7, -2.281087882436601E-4, 3.8658578037323223E-4, 5.2409185656047237E-5, -5.754487394606773E-4, -1.2109014634628852E-4, 7.998633615368073E-5, -7.238929662723885E-5, -1.0736927485933007E-4, 6.86909511654274E-5, 4.896297162612108E-5, -3.004045699599782E-4, 5.8351330755177476E-5, -4.4854846944895313E-4, -2.9324258813092133E-4, 2.9226133030185405E-4, -1.1062183646250851E-4, -6.945484726463543E-4, -2.9176746761405795E-4, -4.3966487452134476E-4, -2.8628069504938163E-4, -2.7048034248184093E-4, 1.1661025996625388E-4, -1.4671765650323154E-4, -4.1704592275336114E-4, -2.371590603077831E-4, 7.742861886695366E-5, -5.587642912810925E-4, -3.9795004925081E-4, -8.814974379197324E-5, 7.396707542698592E-5, 1.5050533015877345E-4, -5.808097598503765E-5, -5.80228622046614E-4, 3.702555949532791E-5, -4.371179800059564E-4, -5.303864236400649E-4, 1.892231237279268E-4, 7.665766038563387E-5, -3.397704931033107E-4, -2.4667218484113667E-4, -1.9284929660263485E-5, -3.032653057142724E-4, 7.819632398016615E-5, -5.172233601295713E-5, -2.5154356171854127E-4, 3.217402779091779E-4, -1.312999796009176E-4, 1.7637771573740123E-4, 2.2168101726111945E-5, -3.1183286687503424E-5, 4.6223926069007795E-5, 1.802637723237964E-4, 3.170084641869661E-4, -1.8924484276680833E-4, -6.946097245888214E-5, 1.1461008485850696E-4, -4.0124588302829133E-4, 7.857072998801834E-5, -2.7324082329048397E-4, -4.1562595212507895E-4, 2.4534202124251367E-4, -8.144624924902777E-5, 6.51632210722416E-5, 3.4389556960025204E-4, -1.761492399835747E-4, 1.1380409311046992E-4, -1.9442410657318578E-4, -4.476621941029252E-4, -2.4750030703671284E-4, 1.1852475520234939E-4, -5.435593818434734E-4, -1.7553337110044103E-4, -1.1124504617850348E-4, -1.9123973293315116E-4, 5.6662982858282683E-5, -1.730531575694268E-4, -8.677986781294714E-5, -2.4326570783304911E-4, 1.789772741660233E-4, 5.670146037351294E-4, -3.713736642547888E-4, 1.2426724156265978E-4, -3.75138672214039E-4, -1.919985513625927E-4, -1.519174311737432E-4, 3.5967806685688114E-5, -2.3051920335899739E-4, 6.925272488579181E-5, 2.3699600710714469E-4, 2.8506025805555523E-4, -4.900881488469072E-4, -6.631807939043097E-5, 4.0315410320179165E-4, -3.888206060628052E-4, 3.0668315554599533E-4, 1.516487479841578E-4, -6.906191398032156E-4, 3.016984300591632E-4, -4.5254932632039985E-4, -3.5160547450370776E-4, -1.0741771557082227E-4, 1.1027304268079502E-4, -6.986027954012375E-4, 4.1755549841171246E-5, -1.2575343619913382E-4, -4.5216318999454137E-4, 1.5476582086842213E-5, -2.5798195672792004E-4, 8.804897924584047E-5, -2.727839908642796E-4, -1.9358308953504842E-5, 3.561535530436805E-4, -1.5637886004566645E-4, 1.224038636873826E-4, 3.2012388394928265E-4, 3.3161805329624816E-4, -5.784691253247822E-4, -1.6741082471132945E-4, -5.220578497127798E-4, -3.405478433809998E-4, -4.113321357431294E-4, 2.5243924055553765E-4, -4.662542331764701E-4, -4.939498375688027E-4, 9.067471357760771E-5, 1.3438552081224998E-4, -2.1796811456036274E-4, -5.571216504120783E-4, -2.1769598840087933E-4, 3.294632249034975E-4, -2.7696897401212903E-4, 1.4597327239297393E-4, 4.106184657462638E-5, 1.6950535796109818E-4, 2.834874582411346E-4, 2.860527117530709E-5, -5.336394127560225E-4, 7.096134384062896E-5, -3.8702190634403333E-4, -3.6475841878189135E-4, -3.895227967248305E-4, 1.1446882499000432E-4, -5.245258237858635E-4, -2.3861746538226375E-4, 3.8896731051558573E-4, -3.3842363776153973E-4, 4.2583393977286264E-4, -1.1364589773183342E-5, -2.178923087303724E-4, 8.587734561686124E-5, -1.4568285997216114E-4, -1.9468386492182813E-4, 8.242753013408587E-5, -4.955510770308594E-4, -7.636504626695056E-5, 3.473214335988908E-4, 1.5900232813095003E-4, -3.326955777346951E-4, 1.0615840897177281E-4, 3.3995660832712607E-4, -5.320802109548995E-4, 2.2149061017945598E-4, -6.57685499671967E-4, -4.348807940353357E-4, 3.6012617339009726E-4, -5.407520183996922E-6, -7.86838371266684E-5, 2.0713951603956303E-4, -7.072630569209135E-4, -3.7351289566815466E-4, -0.0010771183791819165, -5.214541094335232E-4, -3.257325284796718E-4, 3.503699583654106E-5, -2.952908535418219E-4, -3.261800543799613E-4, -2.667797188163929E-4, -5.8557057869901784E-5, -3.3985565866206973E-4, -1.852901360062078E-4, 6.807644372945484E-5, 9.897356592602997E-6, 8.91816753425007E-5, -4.0039833789668117E-4, -1.7640936982793914E-4, 4.5522429030943464E-5, -3.3078676214416095E-4, -7.488133532122191E-4, 2.7738869512226785E-4, -1.5282292548714183E-5, -3.188327369882095E-4, -1.7082236809008733E-4, -2.2219736496775647E-4, -6.26303884259876E-5, -8.579348716922408E-4, -4.154566815416426E-4, 2.825797194660971E-4, -1.1156891975672221E-4, 3.4368019435916555E-4, 6.663494239747501E-5, -8.167947788040667E-4, 4.8398374598097627E-4, -5.711519322985175E-4, -4.119675441203787E-4, -1.497989863770956E-4, 7.90040123949711E-5, -3.730230914612829E-4, -3.09326434268366E-4, -4.359073613795292E-4, 5.202824344008089E-4, 3.2973522897741683E-4, -5.457247034859791E-4, -2.2823944290108246E-4, -2.472196107518218E-4, -2.7829324926909244E-4, -5.153659163792902E-4, 1.3941441025757306E-4, 2.3713990614016516E-4, 1.254710922836945E-4, 1.6061774335608637E-4, -3.8590088462571594E-4, -3.815615330529799E-4, -3.0997576507536763E-4, -1.3724387000644438E-4, -0.0012980608620790871, 9.517590160701368E-5, -3.7650226334284224E-4, -2.5001394818343085E-4, 5.165986854288435E-5, -1.4738146875017347E-4, -5.742522195927131E-4, -9.514995699980804E-4, 4.634495458152788E-4, 1.3403425955448473E-4, -1.0403474162722082E-4, -1.503285286354538E-4, 5.191759329389459E-4, 3.4732614254629705E-4, 4.481199274966334E-4, -1.850377215027905E-4, -9.175974111549066E-5, 2.3411971904778293E-4, -2.729739827908657E-4, -2.657540238657668E-4, -3.992982212836007E-4, 2.068935777158884E-4, -1.0691085526519126E-4, 1.7675976678750645E-4, -5.831937738700406E-5, -2.921838692812118E-4, -2.6270092633652847E-4, 2.04900334399809E-4, 3.614951573787126E-4, -3.4675517561608864E-4, -2.0362192640969437E-4, -1.0213250011505397E-4, 2.4718259083872417E-4, 1.102206393286463E-4, 3.3560476016968656E-4, -8.683189533593787E-5, -1.3725842911153466E-4, -5.68611300005809E-4, -7.062435034042593E-4, -4.273747213757485E-4, -8.853210221847758E-5, 1.298601681454009E-4, -8.254135107436862E-4, -3.6481553462406854E-4, 1.3499020499600323E-4, -3.308740565246887E-5, -6.513067522659998E-4, 3.0694857477190084E-4, 9.910686652920956E-5, 1.8485668149518164E-4, -9.91468793036572E-4, -1.0891998657247895E-5, 2.7975518091384534E-4, 3.585558909843334E-5, 4.032748265219929E-4, 1.036927275604938E-4, -3.0116288143510936E-4, 2.9328560439392245E-4, -8.658578330041698E-4, -3.1927470150517397E-4, -2.4060966081195195E-4, -3.0099637942715955E-4, 5.589101851427575E-5, 2.3914773354824178E-4, -3.0107730613891355E-5, 3.7569920236618344E-4, 3.71849577050133E-4, -2.589608552372837E-4, 3.40824555877121E-4, 1.6223018545442414E-4, -7.616333933987128E-7, -2.477308642150818E-4, 6.762912295399078E-5, -7.539545030652104E-6, 2.3839914874439764E-4, 3.946680914529992E-4, -5.334526667784758E-4, -3.215614042526649E-4, 1.1392266679043908E-4, 5.954570480720613E-5, -2.6809079646941485E-4, -6.029388246762857E-4, -4.3846457836010927E-4, -2.985346225693516E-4, 1.6383243492005628E-4, -9.290877204417115E-5, -4.7301722212781457E-4, -3.4063834686418526E-4, 3.26980131154226E-4, -5.037206711994498E-4, -6.515880961740527E-4, -2.5365009468134177E-4, -5.463349379778966E-4, -8.50764046803989E-5, 2.067018487935426E-4, 1.263611097344144E-4, -7.028387880508784E-4, 2.49745552035742E-4, -7.729187354665556E-4, -1.484697588642933E-4, -2.1490577188519635E-4, 5.748965141506813E-5, -2.453718173516519E-4, -2.25106774404703E-4, 6.313165732909016E-5, -1.1694381145138709E-4, 1.8531751944596324E-4, -2.4361552084469477E-4, -2.4891159732177645E-4, -6.087908799087057E-4, -2.6570229350950743E-4, 7.230666214185486E-5, -3.591460820225956E-5, -1.1946101350387116E-4, -7.992243387831447E-5, 2.2279465400022714E-5, -3.713853551063306E-4, -4.519956758908333E-4, 3.125167430371025E-4, 1.6471310292342427E-4, -8.410948804556636E-4, -1.2166912481419287E-5, -2.689087980134754E-4, -1.5572394985624403E-4, 1.7662238272374948E-4, 6.073606252127845E-5, -3.801101792156574E-4, 1.2715762746974823E-4, -6.668785040681159E-4, -4.40670568475815E-5, -0.0011033257137707461, -3.6277845421011826E-4, -2.0110427270957217E-4, 1.2054074299628226E-4, -3.5465506802002424E-5, -3.297477720070363E-4, -6.509704022022236E-4, 1.453666660011174E-4, -3.2989734195804413E-4, -3.068985477683709E-4, -1.126784825433878E-4, -7.980743007103411E-5, -8.494179744998902E-5, -3.6233726143822384E-4, -2.2593287351741062E-4, 4.052917164781987E-4, -5.410178460320185E-4, -5.793582079395541E-4, 1.8268757851511077E-4, -1.848308838630448E-4, -3.8622214251218253E-4, 3.375021700380606E-4, 1.3904309529194788E-4, -1.7298400746413946E-4, -5.269309211796223E-4, -2.9339028572643246E-4, -1.2763990707792693E-4, -6.603386597320586E-4, 2.583761295810468E-4, -2.973541653706331E-4, -5.11273909730508E-4, -6.050605729228484E-5, -7.091651422961559E-4, -5.555729375147189E-4, -6.649675283049018E-4, 9.497585396750381E-5, -3.590018350986224E-4, 3.9657861078470913E-4, -1.55155249758038E-4, 1.1261188533837562E-4, -2.3114426310888062E-4, 1.8012766439585396E-4, 8.653348051582729E-5, -5.365032060438148E-4, -1.4411798005552112E-4, 1.048817351194945E-4, 9.23584802238171E-5, -1.393928554990765E-4, 4.1958060951416073E-4, 2.3189498911632537E-4, -3.177620781333396E-4, 8.778690530071489E-5, -9.83755366424297E-4, -2.683195481147675E-4, -3.526821233206833E-4, -3.240180316253587E-5, -9.644388173644291E-4, -4.159211245973395E-4, 1.2211181334782772E-4, -2.4323435068373525E-4, -3.438633851973874E-4, -7.646119599648742E-4, -2.298909243909561E-5, 2.1139979923893401E-4, -3.235166415134019E-4, -2.7566882476727793E-4, 4.4107764424213697E-5, 2.270507720683046E-4, 1.2588235442835836E-4, 1.121914167091561E-4, -5.123092258448902E-4, 1.7376848420629808E-4, -7.458962682413239E-4, -3.457577659954678E-4, 1.1761978385687081E-5, 6.040967135640806E-5, -1.7238330277468588E-4, -3.543246509023695E-4, -2.4399415310590717E-4, 1.564911866290858E-4, -8.385259370185312E-5, -3.073078151728803E-5, 4.297569843653425E-5, 3.045768020319149E-4, -1.6199912437694984E-4, 2.3821045986484573E-5, 4.948582290082592E-5, -7.129175505168912E-5, -2.445178729174822E-5, -2.5918990649703103E-5, -5.620415172000923E-4, -3.47042135338144E-4, -2.2542767664111862E-4, 2.1494267773480158E-4, -6.224212133609338E-4, -6.31031925341726E-4, -7.186903144278498E-4, -4.85361112656047E-4, -8.45516036383537E-5, 8.635924686185956E-5, -2.1724964191041423E-4, -4.5851452181828295E-4, -5.087329138425331E-4, -5.2464881988711725E-5, -4.569086751258911E-4, -4.24162230422838E-4, 3.860425800017041E-5, -3.998643686311232E-4, 2.3352507631579385E-4, 4.751535029592002E-5, -4.6463305014013363E-4, 1.7663568252507592E-4, -4.467324852544149E-4, -4.55949839337249E-4, -4.626174901178525E-4, 5.862345370533373E-5, -3.2276289967914627E-4, 1.235671749576085E-4, 1.123696672135957E-4, -3.6720502337148444E-4, 2.2725595573623505E-4, -2.55057134391193E-5, -2.4477569299984986E-4, 2.7946367002761974E-4, -2.976093410841268E-4, -1.2651089067629322E-4, 8.6278124447948E-5, -1.282531837353871E-4, -9.80371640199822E-5, 7.565807843626089E-5, 9.648408608042608E-5, -2.9996295775290126E-4, 2.4163473069715275E-4, 2.527489477459808E-4, -6.772354939098346E-4, 4.3679385839095184E-5, -3.1737934746122083E-4, -4.998878472319067E-4, 1.2756953098582054E-4, -3.0879008251486395E-4, -2.249945638977259E-5, -3.9765226875038324E-5, -4.141488121682292E-4, 1.9927724074285998E-4, -8.66360764141257E-4, -1.59854463798388E-4, -3.9620997780039373E-4, -7.024769154271657E-5, -5.52328266707881E-4, 9.013914535660646E-5, -1.3371652006929138E-4, -3.556572404433717E-5, -3.672870025637565E-4, -2.1118480030032027E-4, 1.3795854583794123E-5, -1.808028782129226E-5, 1.4073407278645557E-4, 3.761708569255418E-4, -3.4384101102210245E-4, -1.5397851112605628E-4, -5.747458371047527E-4, -3.4285365947148237E-5, -3.208972302219612E-4, -1.9172391529535923E-5, -5.16170269967055E-4, -6.493100296740676E-4, -1.2385280243335406E-4, -3.954764106995792E-5, -5.986072107575981E-4, -2.1222033725023343E-4, -3.678298505882787E-5, -1.509159873104521E-4, 3.2669798066979795E-4, 4.098705269970303E-5, -6.285531346962131E-4, 2.4713028841336104E-5, -9.411806093767543E-4, -6.069572935136891E-4, -8.409539146979296E-5, 1.492669675403725E-4, -5.865615170127171E-4, -7.999861520831004E-5, -1.670722355457869E-4, 1.8337453668391263E-6, -1.9786154731778343E-4, 1.706938326229547E-5, 2.4230001041944139E-4, -2.4430482656071775E-4, 2.635559193224199E-4, 3.307702939909878E-5, -3.660406240239714E-5, -2.0400055656709295E-4, 4.34331114773159E-4, 1.703027532954798E-4, -4.886237872026763E-4, -4.9135162283776786E-5, -8.302311368223157E-4, -2.616053672711087E-4, -4.349410558739253E-4, 1.3883261218682746E-4, -6.148798666334855E-4, -1.8776261842721738E-4, -2.0230442890844453E-4, -7.815523610085496E-5, -9.510397640820098E-5, -2.6916514886726013E-4, -2.1655492658882627E-4, 1.7179981799620553E-4, -6.547223623448983E-4, -8.903407023116809E-5, 2.382885565017439E-4, -4.082643158248206E-4, 2.0129388063781144E-4, 3.0320588596971477E-5, -4.5791434477703377E-4, 4.4854888072439364E-4, -6.803469493716287E-4, -2.7108022616751466E-4, -7.495494716845989E-10, 1.3574457628858042E-4, -1.3870972835631427E-4, -1.8521852337849247E-5, 3.884881409934042E-5, 9.255753941804509E-5, 7.993486904749535E-5, 1.8510851297896813E-4, -1.232539824814436E-4, 4.62422819637548E-4, -1.3378373993926703E-4, -1.8610446031235256E-4, 3.79668030143811E-5, -3.329348548421194E-4, 3.112772529968832E-4, 1.3598033457210598E-4, -7.447445682833852E-4, -4.2587633098030366E-4, 1.6186986366592588E-5, 1.5333633438084312E-4, -3.9673663783749004E-4, -5.406352592079363E-4, -2.919516389607192E-4, 2.2215055975559382E-5, -5.808579425964876E-5, 5.2713950710509905E-5, -3.031969716646595E-4, -4.438665798013166E-4, -5.549201106750316E-4, -1.4598879629604375E-4, -8.132974063551411E-4, 1.7993450518053486E-4, 5.680418441332471E-4, -8.382361886533735E-5, -2.8658764461724916E-5, -3.5623004935836985E-5, -1.5018024024764277E-4, 2.2456161230222462E-4, -3.881207769604454E-4, -2.747342878548045E-4, -2.0177082953301605E-4, -2.155005917312793E-4, -4.0914387240057663E-4, 1.2230115966612191E-5, 4.7412849300948186E-4, -2.7878837716190823E-4, 3.64738706513429E-4, -1.1200850344899609E-4, -2.0114817508314216E-4, 2.116726013118395E-4, -2.848933383070848E-4, 4.969672308166384E-4, 3.49931794252441E-4, 1.433769777099928E-4, 2.0735365398155345E-4, 1.4598586070483277E-4, -1.9675209676549917E-6, 5.37118950185325E-5, -5.562362069318975E-5, 2.846822334430818E-4, -4.714161909173229E-4, 2.03376131513277E-4, -5.977838810448473E-4, -3.400526318073431E-4, 2.23715009049113E-4, -3.8260769447297276E-5, 6.47506731151877E-5, -2.523478757813602E-6, -7.44849176196447E-4, -1.9546298305999808E-5, -5.578517324507518E-4, -4.5144082672259425E-4, -1.2130561164641477E-4, 5.2502206759772775E-5, -8.480955534713676E-6, -4.8267910291533377E-4, -1.3615537498764868E-4, 2.7870464466837625E-4, -3.564598310478196E-4, -1.8219287123311993E-4, 1.8254074317463142E-4, 2.456322451501148E-4, 1.6310094073614977E-4, 6.489863738844961E-5, 6.63924097468715E-5, 1.5748177330000677E-4, 1.9826052450348358E-4, -1.4289757151643358E-4, -1.401552675522025E-4, -3.101551467816414E-5, -2.899756409817652E-4, 6.229891540101318E-4, 4.341218955564992E-4, 2.28030105995561E-4, -2.6456980871275267E-4, -2.88483605271464E-4, 4.472085226620064E-5, 1.1521763368358123E-4, 4.481895609370265E-4, -2.2049892861488818E-4, -5.676661124352396E-4, 3.5463870552572617E-4, -2.477454494667166E-4, -2.0370080986689102E-4, -1.5950733069365375E-4, 1.362822685809646E-4, 1.5096862728779616E-5, 3.8969793015948724E-5, -4.068470529322508E-4, -1.205908572763562E-4, 9.656177281891872E-5, 3.359365416662649E-4, 1.9965567542131322E-4, -1.0638573178710169E-4, 2.3877215073202788E-4, 1.4100204123671314E-4, 1.4724691529762637E-4, 2.627746185409279E-4, 1.391374452719891E-4, -1.9345384663460921E-4, -6.025602469825446E-4, -1.95384499426047E-4, -0.0013619118008727023, -4.5462575810018154E-4, -3.9301664512994815E-4, 2.1632020483406295E-4, -2.775685994638614E-4, -9.879054971501078E-5, -4.583831358335288E-5, -1.1126410892539108E-4, -7.255830863251995E-4, 2.0920075620752985E-4, -9.168844723607666E-5, 3.2983377967978024E-4, -3.5142602240706964E-4, -2.1853891266518139E-4, -1.577263721820396E-4, 1.4285534114398145E-4, -4.984608954502947E-5, -1.5620230744305217E-5, -2.1713423262568923E-4, 2.6090308783384983E-4, -9.553175614860849E-4, -5.9754982446439276E-5, -6.390066354934045E-5, 5.510279645616536E-5, -3.392436050686557E-4, 8.822684996504003E-5, 1.324834011116635E-4, 2.4680738655199694E-4, 9.668119371225147E-5, 2.9708045909184345E-4, -2.524898263841875E-4, -7.491994389832246E-5, -2.610254602460288E-5, -3.8259219271414235E-4, 3.660348737032767E-4, 1.931816979612274E-4, 4.747197200225396E-5, -7.5342335088747926E-6, -2.9082944806254176E-4, -9.226429441567977E-5, 3.1837010566254E-4, -9.418989199051113E-5, -2.7550205081024637E-4, -7.035361151900453E-4, -3.059513493944668E-4, -3.3690759825935634E-4, -1.1946857897614238E-4, -8.022364490278822E-5, -4.757390370134784E-4, -2.7703802756857965E-4, 
      };

      System.out.println("the # of ECI are: " + eci.length);

      ce.setECI(eci, activeGroups, false);
      return ce;
	  }
	  
	  
	    /**
	     * Creates a cluster expansion object from a PRIM file.  The cluster expansion object
	     * contains the clusters, ECI, and structure data used to fit the expansion.
	     *
	     * @return The created cluster expansion
	     */
	    public static ClusterExpansion buildCE() {

	        /**
	         * The PRIM file defines primitive cell for the cluster expansion, including
	         * which atoms are allowed at each site.
	         */
	        PRIM prim = new PRIM(PRIM_FILE);

	        /**
	         * By default, it is assumed that we are looking at a system that is periodic in three
	         * dimensions.  If not, you can set the non-periodic dimensions here.  By convention, you
	         * should probably make the third lattice vector non-periodic.  (The third vector has and
	         * index of 2, because counting starts at 1)
	         */
	        prim.setVectorPeriodicity(2, false);

	        /**
	         * The cluster expansion object contains all of the information used to evaluate a
	         * cluster expansion energy model, including all of the clusters and their ECI.  For now,
	         * the ECI values have not been set.
	         */
	        ClusterExpansion ce = new ClusterExpansion(prim);

	        /**
	         * An initial count of how many cluster groups we have, where a cluster group is an
	         * orbit of symmetrically equivalent clusters.
	         */
	        int numTotalClusters = ce.numGroups();

	        /**
	         *  The first argument is the number of sites between clusters, and
	         *  the second is the maximum distance between sites.
	         */
	        ce.appendClusters(1, 1);
	        int numNewClusters = ce.numGroups() - numTotalClusters;
	        numTotalClusters = ce.numGroups();
	        Status.basic("Added " + numNewClusters + " point clusters");

	        ce.appendClusters(2, 4);
	        numNewClusters = ce.numGroups() - numTotalClusters;
	        numTotalClusters = ce.numGroups();
	        Status.basic("Added " + numNewClusters + " pair clusters");

	        ce.appendClusters(3, 2.8);
	        numNewClusters = ce.numGroups() - numTotalClusters;
	        numTotalClusters = ce.numGroups();
	        Status.basic("Added " + numNewClusters + " triple clusters");

	        /**
	         * Print out information about the generated structures of clusters, "Cu" sites
	         * stand for the atoms involved in the certain cluster.
	         */
	        int numFunctionGroups = 0;
	        for (int clustNum = 0; clustNum < ce.numGroups(); clustNum++) {
	          ClusterGroup group = ce.getClusterGroup(clustNum);
	          Status.detail("clustNum:" + clustNum + " , numSites:"
	              + group.numSitesPerCluster() + ", maxdistance:"
	              + group.getMaxDistance() + ", multiplicity:"
	              + group.getMultiplicity()+", ECI:"
	              + group.getECI(0));
	          POSCAR outfile = new POSCAR(group.getSampleStructure(Species.get("Cu"), Species.get("Pt")));
	          outfile.writeFile(CLUSTER_DIR + "cluster." + clustNum + ".vasp");
	          numFunctionGroups += group.numFunctionGroups();
	        }
	        System.out.println(ce.numGroups() + " cluster groups generated");
	        System.out.println(numFunctionGroups + " function groups generated\n\n");

	        return ce;

	      }
	  
}
