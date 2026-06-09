package HEA_fcc_CoIrNiRhRu_211_GitHub;



import matsci.Species;

import matsci.structure.decorate.function.AppliedDecorationFunction;
import matsci.engine.monte.GroundStateRecorder;
import matsci.engine.monte.IAllowsSnapshot;
import matsci.structure.decorate.function.ce.AbstractAppliedCE;

public class CombinedRecorder extends GroundStateRecorder {

	  protected double[][] m_AverageStateCounts;

	  protected AppliedDecorationFunction m_AppliedFunction;
	
	  protected final AbstractAppliedCE m_AppliedCE;

	  protected double m_AverageSOccupancy = 0;
	  protected double m_AverageAdsEnergy = 0;
	  protected double m_AverageActivity = 0;
	  protected double m_AverageAdsEnergySq = 0;
	  protected double m_AverageActivitySq = 0;

	  protected int[] m_SulfurSites;
	  protected int[] m_SulfurStates;
	  protected int[] m_VacancyStates;

	  protected static double m_RuAdsEnergy = -0.875808835;		//TOF benchmark for normalization
	  protected final double m_NChemPot;

	  protected double[] m_RecorderEads40Bars_AdsE;
	  protected double[] m_RecorderHistogram_AdsE;
	  protected double[] m_RecorderEads40Bars_Acti;
	  protected double[] m_RecorderHistogram_Acti;

	  protected int NumSnapshots;

	  
	  
	  
	  public CombinedRecorder(IAllowsSnapshot system, AbstractAppliedCE appliedCE) {
		super(system);

		
		
	    m_AppliedCE = appliedCE;

	    m_NChemPot = appliedCE.getChemPot(Species.nitrogen);

	    int numSSites = 0;
	    for (int sigmaIndex = 0; sigmaIndex < appliedCE.numSigmaSites(); sigmaIndex++) {
	      if (appliedCE.allowsSpecies(sigmaIndex, Species.nitrogen)) {numSSites++;}
	    }

	    m_SulfurSites = new int[numSSites];
	    m_SulfurStates = new int[numSSites];
	    m_VacancyStates = new int[numSSites];

	    int oSiteIndex = 0;
	    for (int sigmaIndex = 0; sigmaIndex < appliedCE.numSigmaSites(); sigmaIndex++) {
	      if (!appliedCE.allowsSpecies(sigmaIndex, Species.nitrogen)) {continue;}
	      m_SulfurSites[oSiteIndex] = sigmaIndex;
	      m_SulfurStates[oSiteIndex] = appliedCE.getStateForSpecies(sigmaIndex, Species.nitrogen);
	      m_VacancyStates[oSiteIndex] = appliedCE.getStateForSpecies(sigmaIndex, Species.vacancy);
	      oSiteIndex++;
	    }

	    m_RecorderEads40Bars_AdsE = new double[40];
	    m_RecorderHistogram_AdsE = new double[40];
	    m_RecorderEads40Bars_Acti = new double[40];
	    m_RecorderHistogram_Acti = new double[40];
	    NumSnapshots = 0;

	    this.m_AppliedFunction = m_AppliedCE;
		this.m_AverageStateCounts = new double[m_AppliedCE.numSigmaSites()][];
		for (int sigmaIndex = 0; sigmaIndex < this.m_AverageStateCounts.length; sigmaIndex++) {
		  int numStates = (m_AppliedCE.getAllowedSpecies(sigmaIndex)).length;
		  this.m_AverageStateCounts[sigmaIndex] = new double[numStates];
		}
	  }
	  
	  static public double EAds2Activity(double EAds) {	//predict TOF value from volcano plot
		  double[] m_EAds2Activity = new double[]{		//Stores the MKM-fitted volcano curve
				  1.87E-06,2.42E-06,3.14E-06,4.06E-06,5.26E-06,6.81E-06,8.82E-06,1.14E-05,1.48E-05,1.92E-05,2.48E-05,3.21E-05,4.16E-05,5.39E-05,6.97E-05,9.03E-05,0.000116965,0.000151466,0.000196144,0.000253999,0.000328918,0.000425933,0.000551562,0.000714239,0.000924892,0.001197663,0.001550866,0.002008211,0.002600391,0.003367141,0.004359892,0.005645215,0.007309265,0.009463527,0.012252245,0.015862014,0.020534166,0.026580738,0.034405075,0.044528346,0.057623662,0.074559913,0.096457977,0.124762676,0.161334647,0.208567296,0.269535161,0.348181253,0.449552375,0.580092721,0.748007183,0.963706335,1.240344489,1.594459754,2.046719579,2.62276544,3.354134379,4.279211479,5.444134798,6.903533164,8.720930772,10.96860844,13.72668371,17.08118072,21.1209293,25.93327992,31.59884954,38.18579321,45.74436255,54.30266845,63.8645213,74.40993527,85.89839837,98.27445824,111.4747238,125.4351711,140.0977128,155.4152816,171.3550775,187.8999976,205.0485359,222.8135708,241.220475,260.3049187,280.1106447,300.6873912,322.0890565,344.372136,367.5944215,391.8139307,417.0880234,443.472661,471.021767,499.7866543,529.8154879,561.152763,593.8387772,627.9090877,663.3939429,700.317685,738.6981215,778.5458668,819.8636574,862.6456468,906.8766879,952.5316127,999.5745201,1047.958087,1097.622918,1148.496943,1200.494895,1253.517872,1307.45301,1362.173285,1417.537464,1473.390222,1529.562441,1585.871714,1642.12305,1698.109808,1753.614855,1808.411941,1862.267301,1914.941469,1966.191278,2015.772034,2063.439838,2108.95401,2152.079593,2192.589876,2230.26892,2264.913999,2296.337952,2324.371365,2348.864562,2369.689344,2386.740461,2399.936772,2409.222079,2414.565606,2415.962142,2413.431826,2407.019586,2396.794268,2382.847454,2365.29202,2344.260456,2319.902991,2292.385569,2261.887706,2228.60029,2192.723342,2154.463797,2114.033328,2071.646248,2027.517527,1981.86093,1934.887311,1886.803065,1837.808754,1788.097912,1737.856023,1687.259681,1636.475913,1585.661674,1534.963485,1484.517215,1434.447998,1384.870257,1335.887836,1287.594214,1240.072806,1193.397314,1147.632139,1102.832826,1059.046546,1016.312594,974.6628969,934.1225353,894.7102548,856.4389766,819.3162949,783.3449601,748.5233441,714.8458852,682.3035118,650.8840423,620.5725622,591.3517762,563.2023364,536.1031473,510.0316464,484.9640627,460.8756534,437.7409197,415.5338022,394.227858,373.7964199,354.212739,335.4501117,317.4819924,300.2820925,283.8244671,268.0835907,253.0344218,238.6524578,224.9137814,211.7950987,199.2737696,187.3278312,175.9360159,165.077762,154.7332206,144.8832564,135.5094444,126.5940624,118.1200796,110.0711406,102.4315473,95.18623652,88.3207551,81.82123188,75.67434725,69.86730031,64.38777389,59.22389772,54.36421011,49.79761853,45.51335941,41.50095768,37.75018662,34.25102847,30.99363661,27.96829984,25.16540947,22.57543007,20.18887456,17.99628427,15.98821478,14.15522802,12.48789116,10.97678261,9.612505169,8.385706399,7.287105677,6.307527453,5.43793972,4.669496541,3.993583201,3.401862351,2.886319406,2.439305397,2.05357558,1.722322236,1.439200418,1.19834572,0.994383573,0.822430025,0.678084348,0.557414226,0.456934545,0.373581029,0.304680048,0.247915946,0.201297135,0.163122084,0.13194614,0.1065499,0.085909686,0.069170466,0.055621389,0.044674007,0.035843107,0.028730041,0.023008347,0.018411479,0.014722392,0.011764787,0.009395774,0.007499775,0.00598347,0.004771635,0.00380372,0.003031052,0.002414552,0.001922876,0.001530909,0.001218546,0.000969703,0.000771522,0.000613731,0.000488131,0.000388177,0.000308647,0.000245381,0.000195061,0.000155044,0.000123224,9.79E-05,7.78E-05,6.18E-05
		  };
		  double Scale = (EAds + 2) * 100;
		  int Scaleint = (int)Scale;
		  double Scaleres = Scale - Scaleint;
		  double Activity = m_EAds2Activity[Scaleint] * (1 - Scaleres) + m_EAds2Activity[Scaleint+1] * Scaleres;
		  return Activity;
	  }
	  /* (non-Javadoc)
	   * @see matsci.engine.monte.IRecorder#recordState(double, double)
	   */
	  @Override
	public void recordEndState(double value, double weight) {
		  
		    double oldTotalWeight = m_TotalWeight;
		    double newTotalWeight = oldTotalWeight + weight;
		    super.recordEndState(value, weight);
		    double correctionRatio = oldTotalWeight / newTotalWeight;
		    double incrementRatio = (1.0D - correctionRatio);

		    int numCo = m_AppliedCE.countSpeciesOccurences(Species.cobalt);
		    int numIr = m_AppliedCE.countSpeciesOccurences(Species.iridium);
		    int numNi = m_AppliedCE.countSpeciesOccurences(Species.nickel);
		    int numRh = m_AppliedCE.countSpeciesOccurences(Species.rhodium);
		    int numRu = m_AppliedCE.countSpeciesOccurences(Species.ruthenium);
		    double latCo = 3.556;
		    double latIr = 3.889;
		    double latNi = 3.556;
		    double latRh = 3.859;
		    double latRu = 3.830;
		    double latavg = (latCo*numCo + latIr*numIr + latNi*numNi + latRh*numRh + latRu*numRu)/(double)(numCo + numIr + numNi + numRh + numRu);
		    double AdsShift = -3.8989*(latavg/3.738 - 1);

		    double SOccupancy = 0;
		    double averageAdsEnergy = 0;
		    double averageActivity = 0;
		    double averageAdsEnergySq = 0;
		    double averageActivitySq = 0;
		    double ActivityRu = EAds2Activity(m_RuAdsEnergy);

		    for (int sigmaIndex = 0; sigmaIndex < this.m_AverageStateCounts.length; sigmaIndex++) {
		        double[] stateCounts = this.m_AverageStateCounts[sigmaIndex];
		        for (int state = 0; state < stateCounts.length; state++)
		          stateCounts[state] = stateCounts[state] * correctionRatio; 
		    } 
		    for (int sigmaIndex = 0; sigmaIndex < this.m_AppliedFunction.numSigmaSites(); sigmaIndex++) {
		        int state = this.m_AppliedFunction.getQuickSigmaState(sigmaIndex);
		        this.m_AverageStateCounts[sigmaIndex][state] = this.m_AverageStateCounts[sigmaIndex][state] + incrementRatio;
		    }
		    
		    for (int x = 0; x < m_RecorderEads40Bars_AdsE.length; x++) {
		        m_RecorderEads40Bars_AdsE[x] = m_RecorderEads40Bars_AdsE[x] * correctionRatio;
		    }
		    for (int x = 0; x < m_RecorderEads40Bars_Acti.length; x++) {
		        m_RecorderEads40Bars_Acti[x] = m_RecorderEads40Bars_Acti[x] * correctionRatio;
		    }

		  	int structnum = 0;

		  	double[] AdsEnergy = new double[3000];
		  	double[] Activity = new double[3000];

		  	for (int nSiteNum = 0; nSiteNum < m_SulfurSites.length/2; nSiteNum++) { //go through nitrogen sites on one side
		  	 	int sigmaIndex = m_SulfurSites[nSiteNum];
		  	 	int nitrogenState = m_SulfurStates[nSiteNum];

		  	 	double delta =0;
		  	 	if (m_AppliedCE.getQuickSigmaState(sigmaIndex) == nitrogenState) { // occupied by nitrogen;
		  	 	    delta = -m_AppliedCE.getDelta(sigmaIndex, m_VacancyStates[nSiteNum]) - m_NChemPot;
		  	 	}
		  	 	else { // occupied by vacancy: Vac--> nitrogen
		  	 		delta = m_AppliedCE.getDelta(sigmaIndex, nitrogenState) - m_NChemPot;

		  	 	}
		  	 	AdsEnergy[nSiteNum]=delta + AdsShift;
		  	 	Activity[nSiteNum]=EAds2Activity(delta + AdsShift) / ActivityRu;

		  	 	double TempAdsE = delta + AdsShift;
		  	 	double TempActi = EAds2Activity(delta + AdsShift) / ActivityRu;
		  	 	averageAdsEnergy += TempAdsE;
		  	 	averageActivity += TempActi;
		  	 	averageAdsEnergySq += TempAdsE * TempAdsE;
		  	 	averageActivitySq += TempActi * TempActi;
		  	 	SOccupancy += 1;

		  	 	structnum++;
		  	}

		  	averageAdsEnergy /= SOccupancy;
		  	averageActivity /= SOccupancy;
		  	averageAdsEnergySq /= SOccupancy;
		  	averageActivitySq /= SOccupancy;
		  	SOccupancy /= m_SulfurSites.length;

		  	int[] NumDistribution_AdsE = new int[40];//histogram of 40 columns
		  	int[] NumDistribution_Acti = new int[40];//histogram of 40 columns

		  	for(int i = 0; i < structnum; i++) {
		  	 	if(Activity[i] < 0.05){
		  	 		NumDistribution_Acti[0]++;
		  	    }
		  	 	else if(Activity[i] < 0.1){
		  	 		NumDistribution_Acti[1]++;
		  	    }
		  	 	else if(Activity[i] < 0.15){
		  	 		NumDistribution_Acti[2]++;
		  	    }
		  	 	else if(Activity[i] < 0.2){
		  	    	NumDistribution_Acti[3]++;
		  	    }
		  	 	else if(Activity[i] < 0.25){
		  	    	NumDistribution_Acti[4]++;
		  	    }
		  	 	else if(Activity[i] < 0.3){
		  	    	NumDistribution_Acti[5]++;
		  	    }
		  	 	else if(Activity[i] < 0.35){
		  	    	NumDistribution_Acti[6]++;
		  	    }
		  	 	else if(Activity[i] < 0.4){
		  	    	NumDistribution_Acti[7]++;
		  	    }
		  	 	else if(Activity[i] < 0.45){
		  	    	NumDistribution_Acti[8]++;
		  	    }
		  	 	else if(Activity[i] < 0.5){
		  	    	NumDistribution_Acti[9]++;
		  	    }
		  	 	else if(Activity[i] < 0.55){
		  	    	NumDistribution_Acti[10]++;
		  	    }
		  	 	else if(Activity[i] < 0.6){
		  	    	NumDistribution_Acti[11]++;
		  	    }
		  	 	else if(Activity[i] < 0.65){
		  	    	NumDistribution_Acti[12]++;
		  	    }
		  	 	else if(Activity[i] < 0.7){
		  	    	NumDistribution_Acti[13]++;
		  	    }
		  	 	else if(Activity[i] < 0.75){
		  	    	NumDistribution_Acti[14]++;
		  	    }
		  	 	else if(Activity[i] < 0.8){
		  	    	NumDistribution_Acti[15]++;
		  	    }
		  	 	else if(Activity[i] < 0.85){
		  	    	NumDistribution_Acti[16]++;
		  	    }
		  	 	else if(Activity[i] < 0.9){
		  	    	NumDistribution_Acti[17]++;
		  	    }
		  	 	else if(Activity[i] < 0.95){
		  	    	NumDistribution_Acti[18]++;
		  	    }
		  	 	else if(Activity[i] < 1.0){
		  	    	NumDistribution_Acti[19]++;
		  	    }
		  	 	else if(Activity[i] < 1.05){
		  	    	NumDistribution_Acti[20]++;
		  	    }
		  	 	else if(Activity[i] < 1.1){
		  	    	NumDistribution_Acti[21]++;
		  	    }
		  	 	else if(Activity[i] < 1.15){
		  	    	NumDistribution_Acti[22]++;
		  	    }
		  	 	else if(Activity[i] < 1.2){
		  	    	NumDistribution_Acti[23]++;
		  	    }
		  	 	else if(Activity[i] < 1.25){
		  	    	NumDistribution_Acti[24]++;
		  	    }
		  	    else if(Activity[i] < 1.3){
		  	    	NumDistribution_Acti[25]++;
		  	    }
		  	    else if(Activity[i] < 1.35){
		  	    	NumDistribution_Acti[26]++;
		  	    }
		  	    else if(Activity[i] < 1.4){
		  	    	NumDistribution_Acti[27]++;
		  	    }
		  	    else if(Activity[i] < 1.45){
		  	    	NumDistribution_Acti[28]++;
		  	    }
		  	    else if(Activity[i] < 1.5){
		  	    	NumDistribution_Acti[29]++;
		  	    }
		  	    else if(Activity[i] < 1.55){
		  	    	NumDistribution_Acti[30]++;
		  	    }
		  	    else if(Activity[i] < 1.6){
		  	    	NumDistribution_Acti[31]++;
		  	    }
		  	    else if(Activity[i] < 1.65){
		  	    	NumDistribution_Acti[32]++;
		  	    }
		  	    else if(Activity[i] < 1.7){
		  	    	NumDistribution_Acti[33]++;
		  	    }
		  	    else if(Activity[i] < 1.75){
		  	    	NumDistribution_Acti[34]++;
		  	    }
		  	    else if(Activity[i] < 1.8){
		  	    	NumDistribution_Acti[35]++;
		  	    }
		  	    else if(Activity[i] < 1.85){
		  	    	NumDistribution_Acti[36]++;
		  	    }
		  	    else if(Activity[i] < 1.9){
		  	    	NumDistribution_Acti[37]++;
		  	    }
		  	    else if(Activity[i] < 1.95){
		  	    	NumDistribution_Acti[38]++;
		  	    }
		  	    else {
		  	    	NumDistribution_Acti[39]++;
		  	    }
		  	}

		  	for(int i = 0; i < structnum; i++) {
		  	 	if(AdsEnergy[i] < -1.45){
		  	 		NumDistribution_AdsE[0]++;
		  	    }
		  	 	else if(AdsEnergy[i] < -1.4){
		  	 		NumDistribution_AdsE[1]++;
		  	    }
		  	 	else if(AdsEnergy[i] < -1.35){
		  	 		NumDistribution_AdsE[2]++;
		  	    }
		  	 	else if(AdsEnergy[i] < -1.3){
		  	    	NumDistribution_AdsE[3]++;
		  	    }
		  	 	else if(AdsEnergy[i] < -1.25){
		  	    	NumDistribution_AdsE[4]++;
		  	    }
		  	 	else if(AdsEnergy[i] < -1.2){
		  	    	NumDistribution_AdsE[5]++;
		  	    }
		  	 	else if(AdsEnergy[i] < -1.15){
		  	    	NumDistribution_AdsE[6]++;
		  	    }
		  	 	else if(AdsEnergy[i] < -1.1){
		  	    	NumDistribution_AdsE[7]++;
		  	    }
		  	 	else if(AdsEnergy[i] < -1.05){
		  	    	NumDistribution_AdsE[8]++;
		  	    }
		  	 	else if(AdsEnergy[i] < -1.0){
		  	    	NumDistribution_AdsE[9]++;
		  	    }
		  	 	else if(AdsEnergy[i] < -0.95){
		  	    	NumDistribution_AdsE[10]++;
		  	    }
		  	 	else if(AdsEnergy[i] < -0.9){
		  	    	NumDistribution_AdsE[11]++;
		  	    }
		  	 	else if(AdsEnergy[i] < -0.85){
		  	    	NumDistribution_AdsE[12]++;
		  	    }
		  	 	else if(AdsEnergy[i] < -0.8){
		  	    	NumDistribution_AdsE[13]++;
		  	    }
		  	 	else if(AdsEnergy[i] < -0.75){
		  	    	NumDistribution_AdsE[14]++;
		  	    }
		  	 	else if(AdsEnergy[i] < -0.7){
		  	    	NumDistribution_AdsE[15]++;
		  	    }
		  	 	else if(AdsEnergy[i] < -0.65){
		  	    	NumDistribution_AdsE[16]++;
		  	    }
		  	 	else if(AdsEnergy[i] < -0.6){
		  	    	NumDistribution_AdsE[17]++;
		  	    }
		  	 	else if(AdsEnergy[i] < -0.55){
		  	    	NumDistribution_AdsE[18]++;
		  	    }
		  	 	else if(AdsEnergy[i] < -0.5){
		  	    	NumDistribution_AdsE[19]++;
		  	    }
		  	 	else if(AdsEnergy[i] < -0.45){
		  	    	NumDistribution_AdsE[20]++;
		  	    }
		  	 	else if(AdsEnergy[i] < -0.4){
		  	    	NumDistribution_AdsE[21]++;
		  	    }
		  	 	else if(AdsEnergy[i] < -0.35){
		  	    	NumDistribution_AdsE[22]++;
		  	    }
		  	 	else if(AdsEnergy[i] < -0.3){
		  	    	NumDistribution_AdsE[23]++;
		  	    }
		  	 	else if(AdsEnergy[i] < -0.25){
		  	    	NumDistribution_AdsE[24]++;
		  	    }
		  	    else if(AdsEnergy[i] < -0.2){
		  	    	NumDistribution_AdsE[25]++;
		  	    }
		  	    else if(AdsEnergy[i] < -0.15){
		  	    	NumDistribution_AdsE[26]++;
		  	    }
		  	    else if(AdsEnergy[i] < -0.1){
		  	    	NumDistribution_AdsE[27]++;
		  	    }
		  	    else if(AdsEnergy[i] < -0.05){
		  	    	NumDistribution_AdsE[28]++;
		  	    }
		  	    else if(AdsEnergy[i] < 0.0){
		  	    	NumDistribution_AdsE[29]++;
		  	    }
		  	    else if(AdsEnergy[i] < 0.05){
		  	    	NumDistribution_AdsE[30]++;
		  	    }
		  	    else if(AdsEnergy[i] < 0.1){
		  	    	NumDistribution_AdsE[31]++;
		  	    }
		  	    else if(AdsEnergy[i] < 0.15){
		  	    	NumDistribution_AdsE[32]++;
		  	    }
		  	    else if(AdsEnergy[i] < 0.2){
		  	    	NumDistribution_AdsE[33]++;
		  	    }
		  	    else if(AdsEnergy[i] < 0.25){
		  	    	NumDistribution_AdsE[34]++;
		  	    }
		  	    else if(AdsEnergy[i] < 0.3){
		  	    	NumDistribution_AdsE[35]++;
		  	    }
		  	    else if(AdsEnergy[i] < 0.35){
		  	    	NumDistribution_AdsE[36]++;
		  	    }
		  	    else if(AdsEnergy[i] < 0.4){
		  	    	NumDistribution_AdsE[37]++;
		  	    }
		  	    else if(AdsEnergy[i] < 0.45){
		  	    	NumDistribution_AdsE[38]++;
		  	    }
		  	    else{
		  	    	NumDistribution_AdsE[39]++;
		  	    }
		  	}
		  	


		    for (int varnum = 0; varnum < m_RecorderEads40Bars_AdsE.length; varnum++) {
		    	m_RecorderEads40Bars_AdsE[varnum] = m_RecorderEads40Bars_AdsE[varnum] + incrementRatio * NumDistribution_AdsE[varnum];

		    }
		    for (int varnum = 0; varnum < m_RecorderEads40Bars_Acti.length; varnum++) {
		    	m_RecorderEads40Bars_Acti[varnum] = m_RecorderEads40Bars_Acti[varnum] + incrementRatio * NumDistribution_Acti[varnum];

		    }


		    double AdsEsum=0;
		    double Actisum=0;
		    for(int i=0; i<40; i++) {
		    	AdsEsum += m_RecorderEads40Bars_AdsE[i];
		    	Actisum += m_RecorderEads40Bars_Acti[i];
		    }
		    for(int i=0; i<40; i++) {
		    	m_RecorderHistogram_AdsE[i]=m_RecorderEads40Bars_AdsE[i]/AdsEsum;
		    	m_RecorderHistogram_Acti[i]=m_RecorderEads40Bars_Acti[i]/Actisum;
			}

		    m_AverageSOccupancy =  (correctionRatio * m_AverageSOccupancy) + (incrementRatio * SOccupancy);
		    m_AverageAdsEnergy =  (correctionRatio * m_AverageAdsEnergy) + (incrementRatio * averageAdsEnergy);
		    m_AverageActivity =  (correctionRatio * m_AverageActivity) + (incrementRatio * averageActivity);
		    m_AverageAdsEnergySq =  (correctionRatio * m_AverageAdsEnergySq) + (incrementRatio * averageAdsEnergySq);
		    m_AverageActivitySq =  (correctionRatio * m_AverageActivitySq) + (incrementRatio * averageActivitySq);

  }

	
	  public double getAverageStateCount(int sigmaIndex, int state) {
		    return this.m_AverageStateCounts[sigmaIndex][state];
	  }

	  public double getAverageSOccupancy() {
		    return m_AverageSOccupancy;
	  }

	  public double getAverageAdsEnergy() {
		    return m_AverageAdsEnergy;
	  }
	  
	  public double getAverageAdsEnergySq() {
		    return m_AverageAdsEnergySq;
	  }
	  
	  public double getAverageActivity() {
		    return m_AverageActivity;
	  }
	  
	  public double getAverageActivitySq() {
		    return m_AverageActivitySq;
	  }
	  
	  public double[] getRecorderEads40Bars_AdsE(){
	      return m_RecorderEads40Bars_AdsE;
	  }

	  public double[] getRecorderHistogram_AdsE(){
	      return m_RecorderHistogram_AdsE;
	  }

	  public double[] getRecorderEads40Bars_Acti(){
	      return m_RecorderEads40Bars_Acti;
	  }

	  public double[] getRecorderHistogram_Acti(){
	      return m_RecorderHistogram_Acti;
	  }

	  

}
