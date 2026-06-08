package methods;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.NumberFormat;

import matsci.Species;
import matsci.io.app.log.Status;
import matsci.io.vasp.POSCAR;
import matsci.location.symmetry.operations.SpaceGroup;
import matsci.location.symmetry.operations.SymmetryOperation;
import matsci.structure.Structure;
import matsci.structure.decorate.DecorationTemplate;
import matsci.structure.decorate.function.ce.ClusterExpansion;
import matsci.structure.decorate.function.ce.FastAppliedCE;
import matsci.structure.superstructure.SuperStructure;
//import reference_class.AdsorptionEnergyRecorder;
import surface.CoIrNiRhRu_fcc_111_fitCE;
import surface.TestRecorder;
import surface.DecorationMSymCanonicalMidFreezeManager;
import surface.DecorationMSymMixMidFreezeManager;

public class runMonteCanonicalTest {

	public static String ROOT_DIR = "ce/CoIrNiRhRu-211-new/";
    public static String GROUND_DIR=ROOT_DIR + "/groundstates/";

    public static void main(String[] args) {
        // TODO Auto-generated method stub

        Status.setLogLevelDetail();
        //Status.setLogLevelBasic();

        //ClusterExpansion ce = CoIrNiRhRu_fcc_111_OH_fitCE.getPreFittedCE_20220216_No_adsorption();
        ClusterExpansion ce = CoIrNiRhRu_fcc_111_fitCE.getPreFittedCE_20250616_num920();
        //ClusterExpansion ce = CoIrNiRhRu_fcc_111_fitCE.buildCE();
    	//CoIrNiRhRu_fcc_111_OH_fitCE.fitCE(ce);

    	String[] abc = new String[args.length];
        for(int i = 0; i < args.length & i < 100; i++){
               abc[i]=args[i];
        }

        //runMonteCanonicalTest(ce, Integer.parseInt(abc[0]), Integer.parseInt(abc[1]), Integer.parseInt(abc[2]), Integer.parseInt(abc[3]),abc[4],abc[5], Integer.parseInt(abc[6]), Integer.parseInt(abc[7]), Integer.parseInt(abc[8]));
        //runMonteCanonicalTest(ce, Integer.parseInt(abc[0]), Integer.parseInt(abc[1]), Integer.parseInt(abc[2]), Integer.parseInt(abc[3]),abc[4],abc[5], Integer.parseInt(abc[6]), Double.MAX_VALUE, Integer.parseInt(abc[7]));
        runMonteCanonicalTest(ce, 5, 0, 0, 16, "920-2000K-504.168.84.84.840-numpass=5000", "Sym-mid_equal-prim=80.vasp", 5000, 2000, 980);
        //runMonteCanonicalTest(ce, 5, 0, 0, 16, "test", "Sym-mid_equal-prim=80.vasp", 50, 3000, 980);
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
		  public static int[] runMonteCanonicalTest(ClusterExpansion ce, int sizex1, int sizex2, int sizey1, int sizey2, String fileDir, String initialStr, int numPass,double startTemp2,int endTemp2) {

		      NumberFormat numberFormat2 = NumberFormat.getNumberInstance();
		      numberFormat2.setMaximumFractionDigits(4);

		      Status.basic("\n\nComposition profiles for CoIrNiRhRu under 3000K to 1000K.");
		      

		      String fileDir2 = GROUND_DIR + fileDir + "/" + "CompositionProfiles-CoIrNiRhRu-N-relax-numPass=" + numPass + "/";
		      String fileDirTemp = fileDir2 + "/" + "intermediateStrs/";

		      
		      File dirFile2  = new  File(fileDir2);
		      if ( ! (dirFile2.exists())  &&   ! (dirFile2.isDirectory())) {
		              boolean  creadok  =  dirFile2.mkdirs();
		      }

		      File dirFileTemp  = new  File(fileDirTemp);
		      if ( ! (dirFileTemp.exists())  &&   ! (dirFileTemp.isDirectory())) {
		              boolean  creadok  =  dirFileTemp.mkdirs();
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

		      int numCo300K = 0;
		      int numIr300K = 0;
		      int numNi300K = 0;
		      int numRh300K = 0;
		      int numRu300K = 0;
		      int numN300K = 0;

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

		      /*
		       * set the superToDirect[][] so as to make the ce is corresponding to the superStructure made of numOfPrims PRIM cells
		       */
		      //superToDirect = ce.getBaseStructure().getDefiningLattice().getCompactSuperToDirect(numOfPrims);
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
		      //concentrationVac = (double)(numVac) / (numN + numVac);
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

		          FileWriter fw = new FileWriter(fileDir2 + "Temp-AdsE-TOF.txt" , false);
		          BufferedWriter bw = new BufferedWriter(fw);
		          bw.write("Temp	NumCo	NumIr	NumNi	NumRh	NumRu	NumN	formE	formCV	avgAdsE	avgTOF	avgAdsE2	avgTOF2\r\n");

		          FileWriter fw1 = new FileWriter(fileDir2 + "SublatticeProfile.txt", false);
		          BufferedWriter bw1 = new BufferedWriter(fw1);
		          
		          FileWriter fw2 = new FileWriter(fileDir2 + "AverageAdsEnergy.txt", false);
		          BufferedWriter bw2 = new BufferedWriter(fw2);
		          
		          FileWriter fw3 = new FileWriter(fileDir2 + "AverageTOF.txt", false);
		          BufferedWriter bw3 = new BufferedWriter(fw3);
		          
		          FileWriter fw4 = new FileWriter(fileDir2 + "SublatticeProfile-1st.txt", false);
		          BufferedWriter bw4 = new BufferedWriter(fw4);

		          FileWriter fw5 = new FileWriter(fileDir2 + "SublatticeProfile-2nd.txt", false);
		          BufferedWriter bw5 = new BufferedWriter(fw5);

		          FileWriter fw6 = new FileWriter(fileDir2 + "SublatticeProfile-3rd.txt", false);
		          BufferedWriter bw6 = new BufferedWriter(fw6);

		          FileWriter fw7 = new FileWriter(fileDir2 + "SublatticeProfile-4th.txt", false);
		          BufferedWriter bw7 = new BufferedWriter(fw7);

		          FileWriter fw8 = new FileWriter(fileDir2 + "SublatticeProfile-5th.txt", false);
		          BufferedWriter bw8 = new BufferedWriter(fw8);
		          
		          FileWriter fw9 = new FileWriter(fileDir2 + "SublatticeProfile-6th.txt", false);
		          BufferedWriter bw9 = new BufferedWriter(fw9);
		          
		          FileWriter fw10 = new FileWriter(fileDir2 + "SublatticeProfile-7th.txt", false);
		          BufferedWriter bw10 = new BufferedWriter(fw10);
		          
		          FileWriter fw11 = new FileWriter(fileDir2 + "SublatticeProfile-8th.txt", false);
		          BufferedWriter bw11 = new BufferedWriter(fw11);
		          
		          FileWriter fw12 = new FileWriter(fileDir2 + "SublatticeProfile-9th.txt", false);
		          BufferedWriter bw12 = new BufferedWriter(fw12);
		          
		          FileWriter fw13 = new FileWriter(fileDir2 + "SublatticeProfile-10th.txt", false);
		          BufferedWriter bw13 = new BufferedWriter(fw13);

		          FileWriter fw14 = new FileWriter(fileDir2 + "SublatticeProfile-11th.txt", false);
		          BufferedWriter bw14 = new BufferedWriter(fw14);
		          
		          FileWriter fw15 = new FileWriter(fileDir2 + "SublatticeProfile-Adsorption.txt", false);
		          BufferedWriter bw15 = new BufferedWriter(fw15);

		          for (int stepNum = 0; stepNum < numSteps; stepNum++) {

		              double temp = startTemp + stepNum * tempIncrement;
		

		              metropolis.setTemp(kb * temp);
//		                 metropolis.s
		              System.out.println();
		              /*
		               * mix ensemble ensemble, the compositions are fixed.
		               */
		              //AdEnergyMSymGCTargetedAverageCurrentFirstMidFreezeManager manager = new AdEnergyMSymGCTargetedAverageCurrentFirstMidFreezeManager(ce, superToDirect, appliedCE, mirrorOp, target);

		              //DecorationCanonicalFirstMidFreezeManager manager = new DecorationCanonicalFirstMidFreezeManager(appliedCE);

		              //relax first and bottom
		              //DecorationCanonicalMidFreezeDoubleEventManager manager = new DecorationCanonicalMidFreezeDoubleEventManager(appliedCE, mirrorOp, targetPeak, targetCompete);
		              //System.out.println("\nthe manager used here is: DecorationMSymCanonicalManager");

		             // DecorationMSymCanonicalManager manager = new DecorationMSymCanonicalManager(appliedCE, mirrorOp);
		              DecorationMSymCanonicalMidFreezeManager manager = new DecorationMSymCanonicalMidFreezeManager(appliedCE, mirrorOp);
		              System.out.println("the manager used here is: DecorationMSymCanonicalMidFreezeManager, mirror sym and freeze middle layer Canonical");


		              //DecorationMSymMixMidFreezeManager manager = new DecorationMSymMixMidFreezeManager(appliedCE, mirrorOp);
		              //System.out.println("the manager used here is: DecorationMSymMixMidFreezeManager, mirror sym and freeze middle layer Mix");

		              //SiteConcentrationMSymCanonicalCurrentRecorder recorder = new SiteConcentrationMSymCanonicalCurrentRecorder(manager, appliedCE);


		              //SiteConcentrationMSymCanonicalCurrentRecorder recorder = new SiteConcentrationMSymCanonicalCurrentRecorder(manager, appliedCE);

		              TestRecorder recorder = new TestRecorder(manager, appliedCE);
		              //SiteConcentrationCanonicalDoubleCurrentRecorder recorder = new SiteConcentrationCanonicalDoubleCurrentRecorder(manager, appliedCE);

		              //SiteConcentrationRecorder recorder = new SiteConcentrationRecorder(manager, appliedCE);
//		                 GroundStateRecorder recorder = new GroundStateRecorder(manager);

		              int numPassTemp = numPasses;

		/*
		              if(temp == 300){
		                  numPassTemp = 10 * numPasses;
		              }
		              else{
		                  numPassTemp = numPasses;
		              }
		*/
		              // just for symmetrical Monte Carlo simulations
		             // Equilibration
		              metropolis.setNumIterations(numPassTemp * appliedCE.numSigmaSites()/10*9);
		              metropolis.runBasic(manager);

		              // Run for real
		              metropolis.setNumIterations(numPassTemp * appliedCE.numSigmaSites()/10);
		              //if (stepNum+1 == numSteps) {
		            	  metropolis.run(manager, recorder);
		              //}else metropolis.runBasic(manager);

		              System.out.println("Temperature: " + metropolis.getTemp());
		              System.out.println("Trigger ratio: " + recorder.getTriggerRatio());
		              System.out.println("Average Energy: " + recorder.getAverageValue());
		              System.out.println("Average Adsorption Energy: " + recorder.getAverageAdsEnergy());
		              System.out.println("Average TOF: " + recorder.getAverageActivity());

		     /*
		              metropolis.setNumIterations(numPasses * appliedCE.numSigmaSites()); // whatever you want
		              GroundStateRecorder recorder = metropolis.runGroundState(manager);

		              double minValue = recorder.getMinValue();

		              manager.setState(recorder.getGroundState());
		     */


		              /*
		               * print out the ground state for every temperature step
		               *
		               */
		              SuperStructure groundState3 = appliedCE.getSuperStructure();
		              Structure primGroundState3=groundState3.findPrimStructure().getCompactStructure();
		              //Structure primGroundState3=groundState3;
		              int numCoTem = primGroundState3.numDefiningSitesWithSpecies(Species.cobalt);
		              int numIrTem = primGroundState3.numDefiningSitesWithSpecies(Species.iridium);
		              int numNiTem = primGroundState3.numDefiningSitesWithSpecies(Species.nickel);
		              int numRhTem = primGroundState3.numDefiningSitesWithSpecies(Species.rhodium);
		              int numRuTem = primGroundState3.numDefiningSitesWithSpecies(Species.ruthenium);
		              int numNTem = primGroundState3.numDefiningSitesWithSpecies(Species.nitrogen);
		              //int numVacTem = primGroundState3.numDefiningSitesWithSpecies(Species.vacancy);
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
		             double AdsESq = avgAdsE2 - avgAdsE * avgAdsE;
		             double AdsEcV = AdsESq / (kb * temp * temp);
		             AdsEcV /= appliedCE.numSigmaSites();
		             double ActiSq = avgActi2 - avgActi * avgActi;
		             double ActicV = ActiSq / (kb * temp * temp);
		             ActicV /= appliedCE.numSigmaSites();
		             System.out.println("Temp: " + temp + "\tNumPass:" + numPassTemp + "\tConcentration: " + concentrationCo + " " + concentrationIr + " " + concentrationNi + " " + concentrationRh + " " + concentrationRu + " " + concentrationN  + "\tformE: " + calculatedEnergy + "\tAdsE: " + avgAdsE + "\tcv: " + AdsEcV + "\tTOF: " + avgActi + "\tcv: " + ActicV);
		             System.out.println("the num of steps is: " + stepNum);


		             String stringNum = temp + "     " + numCoTem + "     " + numIrTem + "     " + numNiTem + "     " + numRhTem + "     " + numRuTem + "     " + numNTem + "     " + calculatedEnergy +  "     " + FormESq + "     " + avgAdsE + "     " + avgActi + "     " + AdsEcV + "     " + ActicV + "\r\n";
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

		   /*
		    structMetropolis.setNumIterations(numPasses * appliedCE.numSigmaSites()); // whatever you want
		     GroundStateRecorder recorder = structMetropolis.runGroundState(manager);

		     double minValue = recorder.getMinValue();

		     manager.setState(recorder.getGroundState());
		     SuperStructure groundState = appliedCE.getSuperStructure();

		   */

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
		      outfile2.writeFile(fileDir2 + "prim=" + numPrims + "-" + sizex1 + "_" + sizex2 + "_" + sizey1 + "_" + sizey2 + "."+"groundState." + numCoF + "." + numIrF + "." + numNiF + "." + numRhF + "." + numRuF + "-" + numNF + ".vasp");
		      outfile2.writeVICSOutFile(fileDir2 + "prim=" + numPrims + "-" + sizex1 + "_" + sizex2 + "_" + sizey1 + "_" + sizey2 + "."+"groundState." + numCoF + "." + numIrF + "." + numNiF + "." + numRhF + "." + numRuF + "-" + numNF + ".out");

		      String fileDir3 = GROUND_DIR + fileDir + "/";outfile2.writeFile(fileDir3 + numberFormat2.format(chemPotCo) + "." + numberFormat2.format(chemPotIr) + "." + numberFormat2.format(chemPotNi) + "." + numberFormat2.format(chemPotRh) + "." + numberFormat2.format(chemPotRu) + "-" + sizex1 + "_" + sizex2 + "_" + sizey1 + "_" + sizey2 + "."+"groundState." +numCoF + "." + numIrF + "." + numNiF + "." + numRhF + "." + numRuF + "-" + numNF + ".vasp");
		      outfile2.writeVICSOutFile(fileDir3 + numberFormat2.format(chemPotCo) + "." + numberFormat2.format(chemPotIr) + "." + numberFormat2.format(chemPotNi) + "." + numberFormat2.format(chemPotRh) + "." + numberFormat2.format(chemPotRu) + "-" + sizex1 + "_" + sizex2 + "_" + sizey1 + "_" + sizey2 + "."+"groundState." +numCoF + "." + numIrF + "." + numNiF + "." + numRhF + "." + numRuF + "-" + numNF + ".out");

		      //numberFormat.format(chemPotential)

		      int[] numCoIrNiRhRu =  new int[12];
		       numCoIrNiRhRu[0] = numCoF;
		       numCoIrNiRhRu[1] = numIrF;
		       numCoIrNiRhRu[2] = numNiF;
		       numCoIrNiRhRu[3] = numRhF;
		       numCoIrNiRhRu[4] = numRuF;
		       numCoIrNiRhRu[5] = numNF;

		       numCoIrNiRhRu[6] = numCo300K;
		       numCoIrNiRhRu[7] = numIr300K;
		       numCoIrNiRhRu[8] = numNi300K;
		       numCoIrNiRhRu[9] = numRh300K;
		       numCoIrNiRhRu[10] = numRu300K;
		       numCoIrNiRhRu[11] = numN300K;

		       return numCoIrNiRhRu;
		  } // end of runMonte method


}
