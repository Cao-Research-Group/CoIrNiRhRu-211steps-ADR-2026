/**
 * 
 */
package HEA_fcc_CoIrNiRhRu_211_GitHub;

/**
 * @author Liang
 * 
 * January 22 2013,
 *  based on Tim's DecorationCrossCanonicalManager class
 *  
 *  attention:  another canonical ensemble manager class: DecorationCanonicalManager class,
 *  that forces swaps to stay within sublattices.
 *
 */

import java.util.BitSet;
import java.util.Random;


import matsci.Species;
import matsci.engine.monte.metropolis.IAllowsMetropolis;
import matsci.engine.monte.metropolis.IMetropolisEvent;
import matsci.io.vasp.POSCAR;
import matsci.location.Coordinates;
import matsci.location.basis.AbstractBasis;
import matsci.location.symmetry.operations.SpaceGroup;
import matsci.location.symmetry.operations.SymmetryOperation;
import matsci.structure.Structure;
import matsci.structure.StructureBuilder;
import matsci.structure.decorate.function.AppliedDecorationFunction;
import matsci.util.arrays.ArrayUtils;

public class DecorationMSymCanonicalManager implements IAllowsMetropolis {

  protected AppliedDecorationFunction m_AppliedFunction;
  protected Random m_Generator;
  protected Event m_Event;
  protected DoubleSymEvents m_DoubleEvent;
  //protected SymmetryOperation m_SymOp;
  protected int[] m_SigmaIndexMap;

  
  protected Species[] m_SigmaSpecies;
  protected boolean[][] m_AllowedIndicesBySpecies;
  protected int[] m_SpeciesIndexBySigmaIndex;
  
  
  /*
  public DecorationMSymCanonicalManager(AppliedDecorationFunction ah) {
    m_AppliedFunction = ah;
    m_Generator = new Random();
    m_SigmaSpecies =new Species[0];
    m_AllowedIndicesBySpecies = new boolean[0][];
    m_Event = new Event();
    m_DoubleEvent = new DoubleSymEvents();
    
    this.refreshFromStructure();
  }
 */ 
  
  public DecorationMSymCanonicalManager(AppliedDecorationFunction ah, SymmetryOperation symOp) {
      m_AppliedFunction = ah;
      m_Generator = new Random();
      m_SigmaSpecies =new Species[0];
      m_AllowedIndicesBySpecies = new boolean[0][];
      m_Event = new Event();
      m_DoubleEvent = new DoubleSymEvents();
      //m_SymOp = symOp;
      m_SigmaIndexMap = new int[m_AppliedFunction.numSigmaSites()];
      for (int sigmaIndex = 0; sigmaIndex < m_SigmaIndexMap.length; sigmaIndex++) {
        Coordinates siteCoords = m_AppliedFunction.getSigmaSite(sigmaIndex).getCoords();
        Coordinates oppedCoords = symOp.operate(siteCoords);
        Structure.Site oppedSite = m_AppliedFunction.getSuperStructure().getSite(oppedCoords);
        m_SigmaIndexMap[sigmaIndex] = m_AppliedFunction.getSigmaIndex(oppedSite.getIndex());
      }
      this.refreshFromStructure();
      
      System.out.println("\n\n\n****************************************************************************\n\n");
      System.out.println("finish initializing manager: DecorationMSymCanonicalManager");
      System.out.println("\n\n\n****************************************************************************\n\n");
    }
  
  public void refreshFromStructure() {
      
    //System.out.println("m_AppliedFunction.numSigmaSites(): " + m_AppliedFunction.numSigmaSites());  
    int numSigmaSites = m_AppliedFunction.numSigmaSites();
    
    m_SpeciesIndexBySigmaIndex = new int[numSigmaSites];
    for (int sigmaIndex = 0; sigmaIndex < m_SpeciesIndexBySigmaIndex.length; sigmaIndex++) {
      Species[] allowedSpecies = m_AppliedFunction.getAllowedSpecies(sigmaIndex);
      Species siteSpecies = allowedSpecies[m_AppliedFunction.getQuickSigmaState(sigmaIndex)];
      
      //System.out.println(sigmaIndex + ": spec=" + siteSpecies.toString());
      
      int speciesIndex = ArrayUtils.findIndex(m_SigmaSpecies, siteSpecies);
      if (speciesIndex < 0) {
        speciesIndex = m_SigmaSpecies.length; 
        m_SigmaSpecies = (Species[]) ArrayUtils.appendElement(m_SigmaSpecies, siteSpecies);
        m_AllowedIndicesBySpecies = ArrayUtils.appendElement(m_AllowedIndicesBySpecies, new boolean[numSigmaSites]);
        for (int allowedSpecIndex = 0; allowedSpecIndex < numSigmaSites; allowedSpecIndex++) {
          if (m_AppliedFunction.allowsSpecies(allowedSpecIndex, siteSpecies)) {
            m_AllowedIndicesBySpecies[speciesIndex][allowedSpecIndex] = true;
          }
        }
      }
      m_SpeciesIndexBySigmaIndex[sigmaIndex] = speciesIndex;
    }

  }
  
  public AppliedDecorationFunction getAppliedHamiltonian() {
    return m_AppliedFunction;
  }

  /* (non-Javadoc)
   * @see matsci.engine.monte.metropolis.IAllowsMetropolis#getValue()
   */
  public double getValue() {
    return m_AppliedFunction.getValue();
  }

  /* (non-Javadoc)
   * @see matsci.engine.monte.metropolis.IAllowsMetropolis#getRandomEvent()
   */
  public IMetropolisEvent getRandomEvent() {
    
    boolean isSwapAllowed = false;
    boolean isTwoSwapAllowed = false;
    boolean isTwoPointsSym_1 = false;
    boolean isTwoPointsSym_2 = false;
    
    int sigmaIndex1 = -1;
    int sigmaIndex2 = -1;
    int mirrorSymSigmaIndex1 = -1;
    int mirrorSymSigmaIndex2 = -1;
    
    int newSpeciesIndex1=-1;
    int newSpeciesIndex2=-1;
    int newSpeciesIndex2_1=-1;
    int newSpeciesIndex2_2=-1;

    // comment triedSwaps out because:  it's a sanity check for situations in which there might not be any allowed swaps, which is rare.
    //BitSet triedSwaps = new BitSet();
    int numPossibleSwaps = m_SpeciesIndexBySigmaIndex.length * m_SpeciesIndexBySigmaIndex.length;
     
    while (!isSwapAllowed) {
/*        
      if (triedSwaps.cardinality() == numPossibleSwaps) {
        throw new RuntimeException("All swaps have been tried and none are allowed!");
      }
*/
      sigmaIndex1 = m_Generator.nextInt(m_SpeciesIndexBySigmaIndex.length);
      sigmaIndex2 = m_Generator.nextInt(m_SpeciesIndexBySigmaIndex.length);
      int swapIndex = sigmaIndex2 * m_SpeciesIndexBySigmaIndex.length + sigmaIndex1;
      
      //triedSwaps.set(swapIndex);
      
      newSpeciesIndex1 = m_SpeciesIndexBySigmaIndex[sigmaIndex2];
      newSpeciesIndex2 = m_SpeciesIndexBySigmaIndex[sigmaIndex1];
      
      isSwapAllowed = (newSpeciesIndex1 != newSpeciesIndex2); // Don't swap the same species
      isSwapAllowed &= m_AllowedIndicesBySpecies[newSpeciesIndex1][sigmaIndex1] && m_AllowedIndicesBySpecies[newSpeciesIndex2][sigmaIndex2];
 
      
      /*
      Structure objectStr = m_AppliedFunction.getStructure();
      Coordinates oCoords = objectStr.getSiteCoords(sigmaIndex1);
      Coordinates mirrorSymCoords = m_SymOp.operate(oCoords);
      mirrorSymSigmaIndex1 = objectStr.getDefiningSite(mirrorSymCoords).getIndex();
//      System.out.println("print out the mirror symmetrical point index of first point of swap" + mirrorSymSigmaIndex1);
      Coordinates oCoords2 = objectStr.getSiteCoords(sigmaIndex2);
      Coordinates mirrorSymCoords2 = m_SymOp.operate(oCoords2);    
      mirrorSymSigmaIndex2 = objectStr.getDefiningSite(mirrorSymCoords2).getIndex(); */
      
      mirrorSymSigmaIndex1 = m_SigmaIndexMap[sigmaIndex1];
      mirrorSymSigmaIndex2 = m_SigmaIndexMap[sigmaIndex2];
      int swapIndex2 = mirrorSymSigmaIndex2 * m_SpeciesIndexBySigmaIndex.length + mirrorSymSigmaIndex1;
      //triedSwaps.set(swapIndex2);
            
      newSpeciesIndex2_1 = m_SpeciesIndexBySigmaIndex[mirrorSymSigmaIndex2];
      newSpeciesIndex2_2 = m_SpeciesIndexBySigmaIndex[mirrorSymSigmaIndex1];
      
      isSwapAllowed &= ( newSpeciesIndex2_1 != newSpeciesIndex2_2);; // Don't swap the same species
      isSwapAllowed &= m_AllowedIndicesBySpecies[newSpeciesIndex2_1][mirrorSymSigmaIndex1] && m_AllowedIndicesBySpecies[newSpeciesIndex2_2][mirrorSymSigmaIndex2];    
      
      // to decide whether the two sym points have the same species.
      isTwoPointsSym_1 = (newSpeciesIndex1 == newSpeciesIndex2_1);
      isTwoPointsSym_2 = (newSpeciesIndex2 == newSpeciesIndex2_2);
      isSwapAllowed &= isTwoPointsSym_1;
      isSwapAllowed &= isTwoPointsSym_2;      
      

      isTwoSwapAllowed = ((sigmaIndex1 != mirrorSymSigmaIndex1) && (sigmaIndex2 != mirrorSymSigmaIndex2) && (sigmaIndex1 != mirrorSymSigmaIndex2) && (sigmaIndex2 != mirrorSymSigmaIndex1) && (sigmaIndex1 != sigmaIndex2) && (mirrorSymSigmaIndex1 != mirrorSymSigmaIndex2));
     // isTwoSwapAllowed = ((sigmaIndex1 != mirrorSymSigmaIndex1) && (sigmaIndex2 != mirrorSymSigmaIndex2));
      

      isSwapAllowed &= isTwoSwapAllowed; 
    }        
    
    //System.out.println("two allowed swaps: " + "(" + sigmaIndex1 + "," + sigmaIndex2 +"), " + "(" +  mirrorSymSigmaIndex1 + "," + mirrorSymSigmaIndex2 + "); " + "new species:  (" +  newSpeciesIndex1 + "," + newSpeciesIndex2 +"), " + "(" + newSpeciesIndex2_1 + "," + newSpeciesIndex2_2 + ") ");
    
    m_DoubleEvent.setSwap(sigmaIndex1, sigmaIndex2, mirrorSymSigmaIndex1, mirrorSymSigmaIndex2);
    return m_DoubleEvent;
  }

  /* (non-Javadoc)
   * @see matsci.engine.monte.IAllowsSnapshot#setState(java.lang.Object)
   */
  public void setState(Object snapshot) {
    m_AppliedFunction.setSigmaStates((int[]) snapshot);
    this.refreshFromStructure();
  }

  /* (non-Javadoc)
   * @see matsci.engine.monte.IAllowsSnapshot#getSnapshot(java.lang.Object)
   */
  public Object getSnapshot(Object template) {
    return m_AppliedFunction.getQuickSigmaStates((int[]) template);
  }
  
  public class Event implements IMetropolisEvent {

    private int[] m_SigmaIndicesToSwap = new int[2];
    private int[] m_NewSiteStates = new int[2];
    private int[] m_OldSiteStates = new int[2];
    
    private int m_NewSpeciesIndex1;
    private int m_NewSpeciesIndex2;
    
    public void setSwap(int sigmaIndex1, int sigmaIndex2) {
      
      m_NewSpeciesIndex1 = m_SpeciesIndexBySigmaIndex[sigmaIndex2];
      m_NewSpeciesIndex2 = m_SpeciesIndexBySigmaIndex[sigmaIndex1];
      
      m_SigmaIndicesToSwap[0] = sigmaIndex1;
      m_SigmaIndicesToSwap[1] = sigmaIndex2;
      m_NewSiteStates[0] = m_AppliedFunction.getStateForSpecies(sigmaIndex1, m_SigmaSpecies[m_NewSpeciesIndex1]);
      m_NewSiteStates[1] = m_AppliedFunction.getStateForSpecies(sigmaIndex2, m_SigmaSpecies[m_NewSpeciesIndex2]);
        m_OldSiteStates[0] = m_AppliedFunction.getQuickSigmaState(sigmaIndex1);
        m_OldSiteStates[1] = m_AppliedFunction.getQuickSigmaState(sigmaIndex2);
     
        //System.out.println("print out the site states: (" + m_NewSiteStates[0] + ", " + m_NewSiteStates[1] + ").  (" + m_OldSiteStates[0] + ",  " + m_OldSiteStates[1] + ")");
    }
    
      /* (non-Javadoc)
       * @see matsci.engine.monte.metropolis.IMetropolisEvent#trigger()
       */
      public void trigger() {
        
        m_SpeciesIndexBySigmaIndex[m_SigmaIndicesToSwap[0]] = m_NewSpeciesIndex1;
        m_SpeciesIndexBySigmaIndex[m_SigmaIndicesToSwap[1]] = m_NewSpeciesIndex2;
        
        m_AppliedFunction.setSigmaStates(m_SigmaIndicesToSwap, m_NewSiteStates);
      }
    
      /* (non-Javadoc)
       * @see matsci.engine.monte.metropolis.IMetropolisEvent#getDelta()
       */
      public double getDelta() {
        
        return m_AppliedFunction.getDelta(m_SigmaIndicesToSwap, m_NewSiteStates);

      }
    
      /* (non-Javadoc)
       * @see matsci.engine.monte.metropolis.IMetropolisEvent#triggerGetDelta()
       */
      public double triggerGetDelta() {
        
        m_SpeciesIndexBySigmaIndex[m_SigmaIndicesToSwap[0]] = m_NewSpeciesIndex1;
        m_SpeciesIndexBySigmaIndex[m_SigmaIndicesToSwap[1]] = m_NewSpeciesIndex2;
        
        return m_AppliedFunction.setSigmaStatesGetDelta(m_SigmaIndicesToSwap, m_NewSiteStates);
      }
    
      /* (non-Javadoc)
       * @see matsci.engine.monte.metropolis.IMetropolisEvent#reverse()
       */
      public void reverse() {
        
        m_SpeciesIndexBySigmaIndex[m_SigmaIndicesToSwap[0]] = m_NewSpeciesIndex2;
        m_SpeciesIndexBySigmaIndex[m_SigmaIndicesToSwap[1]] = m_NewSpeciesIndex1;
        
        m_AppliedFunction.setSigmaStates(m_SigmaIndicesToSwap, m_OldSiteStates);
      }
      
      
      public int[] getIndicesToSwap(){
          return this.m_SigmaIndicesToSwap;
      }
  }


  public class DoubleSymEvents implements IMetropolisEvent {

      private Event event1;
      private Event event2;
      
      public DoubleSymEvents() {
          this.event1 = new Event();
          this.event2 = new Event();
      }
      
      
      public void setSwap(int sigmaIndex1_1, int sigmaIndex1_2, int sigmaIndex2_1, int sigmaIndex2_2) {

          this.event1.setSwap(sigmaIndex1_1, sigmaIndex1_2);
          this.event2.setSwap(sigmaIndex2_1, sigmaIndex2_2);

      }
      
        /* (non-Javadoc)
         * @see matsci.engine.monte.metropolis.IMetropolisEvent#trigger()
         */
        public void trigger() {
            this.event1.trigger();
            this.event2.trigger();
        }
      
        /* (non-Javadoc)
         * @see matsci.engine.monte.metropolis.IMetropolisEvent#getDelta()
         */
        public double getDelta() {
          
          double energy0 = m_AppliedFunction.getValue();
          
          this.event1.trigger();
          this.event2.trigger();
          double energy1 = m_AppliedFunction.getValue();
          double deltaAfterTwoSwap = energy1 - energy0;
          this.event2.reverse();
          this.event1.reverse();
//          System.out.println("print out to see whether getValue() is energy: " + energy0 + ", " + energy1);
          
          return  deltaAfterTwoSwap;


        }
      
        /* (non-Javadoc)
         * @see matsci.engine.monte.metropolis.IMetropolisEvent#triggerGetDelta()
         */
        public double triggerGetDelta() {
            double energy0 = m_AppliedFunction.getValue();
//            m_AppliedFunction.get
            this.event1.trigger();
            this.event2.trigger();
            double energy1 = m_AppliedFunction.getValue();
            
 //           System.out.println("print out to see whether getValue() is energy: " + energy0 + ", " + energy1);
            
            return (energy1 - energy0);
        }
      
        /* (non-Javadoc)
         * @see matsci.engine.monte.metropolis.IMetropolisEvent#reverse()
         */
        public void reverse() {
            this.event2.reverse();
            this.event1.reverse();
        }
        
        public IMetropolisEvent getEvent1(){
            return this.event1;
            
        }
        
        public int[] getIndicesToSwap(){
            return this.event1.m_SigmaIndicesToSwap;
        }
    }

  
}




/* edited by Liang to make the Monte Carlo simulation swaps symmetrical to make
 * sure the structures after Monter Carlo is symmetrical(parallel to the surface of 2-dimension system)
 * @see matsci.engine.monte.metropolis.IAllowsMetropolis#getRandomEvent()
 */
/*  
public IMetropolisEvent[] getRandomEvent(SymmetryOperation mirrorOp) {
  
  int sigmaIndex1 = m_Generator.nextInt(m_SpeciesIndexBySigmaIndex.length);
  int sigmaIndex2 = m_Generator.nextInt(m_SpeciesIndexBySigmaIndex.length);

  m_Event[0].setSwap(sigmaIndex1, sigmaIndex2);
  
 
  

  SymmetryOperation op = mirrorOp;
  
  Structure objectStr = this.m_AppliedFunction.getStructure();
  
  Coordinates oCoords = objectStr.getSiteCoords(sigmaIndex1);
  

  Coordinates mirrorSymCoords = op.operate(oCoords);

  
  int mirrorSymSigmaIndex1 = objectStr.getDefiningSite(mirrorSymCoords).getIndex();
  System.out.println("print out the mirror symmetrical point index of first point of swap" + mirrorSymSigmaIndex1);
//  StructureBuilder builder = new StructureBuilder(objectStr);
  
  
  Coordinates oCoords2 = objectStr.getSiteCoords(sigmaIndex2);
  Coordinates mirrorSymCoords2 = op.operate(oCoords2);    
  int mirrorSymSigmaIndex2 = objectStr.getDefiningSite(mirrorSymCoords2).getIndex();
  System.out.println("print out the mirror symmetrical point index of second point of swap: " + mirrorSymSigmaIndex2);    
  
  m_Event[1].setSwap(mirrorSymSigmaIndex1, mirrorSymSigmaIndex2);
  
  
  return m_Event;
}
*/ 
