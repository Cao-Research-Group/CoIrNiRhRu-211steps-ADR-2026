package HEA_fcc_CoIrNiRhRu_211_GitHub;


/**
 * 
 */


/**
 * @author Liang
 * 
 * January 22 2013
 *
 */

import matsci.Species;
import matsci.engine.monte.metropolis.IMetropolisEvent;
import matsci.location.symmetry.operations.SymmetryOperation;
import matsci.structure.decorate.DecorationTemplate;
import matsci.structure.decorate.function.AppliedDecorationFunction;


/**
 * @author Tim Mueller
 *
 * Allows swaps between different sublattices
 * 
 * Might generate a lot of failed (disallowed) swaps.
 * 
 */
public class DecorationMSymCanonicalMidFreezeManager extends DecorationMSymCanonicalManager {
  
  protected int m_MiddleStartLimit ;//excluded
  protected int m_MiddleEndLimit;//included
  
  public DecorationMSymCanonicalMidFreezeManager(AppliedDecorationFunction ah, SymmetryOperation symOp) {
      super(ah, symOp);
      
      this.initialFrozenSites();
      
      System.out.println("\n\n\n****************************************************************************\n\n");
      System.out.println("finish initializing manager: DecorationMSymCanonicalMidFreezeManager");
      System.out.println("\n\n\n****************************************************************************\n\n");
    }
  

  /**
   * assume: # of metal layers: 2(n-2)+1=2n-3
   *         # of adsorbates: 2*1=2
   *         # of total sub lattice: (n-2)+1+1=n
   *         
   *         start allowed layer: 3
   *         end   allowed layer: 3
   */
  public void initialFrozenSites() {
      
      int numPrims = m_AppliedFunction.numPrimCells();  

      DecorationTemplate decorTemplate = this.m_AppliedFunction.getHamiltonian();
      
      int totalNumSites = 0;
              
      for(int sublatticeIndex = 0; sublatticeIndex< decorTemplate.numSublattices(); sublatticeIndex++ ){
          int[] sigmaSiteForSublattice = this.m_AppliedFunction.getSigmaIndicesForSublattice(sublatticeIndex);          
          if (this.m_AppliedFunction.allowsSpecies(sigmaSiteForSublattice[0], Species.nitrogen)) {
              continue; // skip oxygen/vacancy layers
          }
         
          int numSitesForSublattice = sigmaSiteForSublattice.length;
          totalNumSites += numSitesForSublattice;
      }
      
      int totalNumLayers = totalNumSites / numPrims;
      
      this.m_MiddleStartLimit = numPrims * 10 - 1; //excluded    //2025.03.26 by caowang fcc(211) 7 layers fix mid one
      this.m_MiddleEndLimit =m_MiddleStartLimit + numPrims*2; //included
      
      System.out.println("Initializing:  ");
      System.out.println("decorTemplate.numSublattices(), m_AppliedFunction.numPrimCells():  " + decorTemplate.numSublattices() + ", " + m_AppliedFunction.numPrimCells());
      System.out.println("total # of layers:  " + totalNumLayers);
      System.out.println("m_MiddleStartLimit:  " + m_MiddleStartLimit);
      System.out.println("m_MiddleEndLimit:  " + m_MiddleEndLimit);
 
  }
  

  /* (non-Javadoc)
   * @see matsci.engine.monte.metropolis.IAllowsMetropolis#getRandomEvent()
   */
  /**
   * freeze middle layer, for example, Pt3Ni, to mimic the bulk core
   */
  public IMetropolisEvent getRandomEvent() {
    
    boolean NotMiddleLayers = false;
      
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
     
    while (!isSwapAllowed) {

      sigmaIndex1 = m_Generator.nextInt(m_SpeciesIndexBySigmaIndex.length);
      sigmaIndex2 = m_Generator.nextInt(m_SpeciesIndexBySigmaIndex.length);
      
      newSpeciesIndex1 = m_SpeciesIndexBySigmaIndex[sigmaIndex2];
      newSpeciesIndex2 = m_SpeciesIndexBySigmaIndex[sigmaIndex1];
      
      isSwapAllowed = (newSpeciesIndex1 != newSpeciesIndex2); // Don't swap the same species
      isSwapAllowed &= m_AllowedIndicesBySpecies[newSpeciesIndex1][sigmaIndex1] && m_AllowedIndicesBySpecies[newSpeciesIndex2][sigmaIndex2];
      
      mirrorSymSigmaIndex1 = m_SigmaIndexMap[sigmaIndex1];
      mirrorSymSigmaIndex2 = m_SigmaIndexMap[sigmaIndex2];
      
      newSpeciesIndex2_1 = m_SpeciesIndexBySigmaIndex[mirrorSymSigmaIndex2];
      newSpeciesIndex2_2 = m_SpeciesIndexBySigmaIndex[mirrorSymSigmaIndex1];
      
      // to decide whether the two sym points have the same species.
      isTwoPointsSym_1 = newSpeciesIndex1 == newSpeciesIndex2_1;
      isTwoPointsSym_2 = newSpeciesIndex2 == newSpeciesIndex2_2;
      isSwapAllowed &= isTwoPointsSym_1;
      isSwapAllowed &= isTwoPointsSym_2;      
      
      isSwapAllowed &= ( newSpeciesIndex2_1 != newSpeciesIndex2_2);; // Don't swap the same species
      isSwapAllowed &= m_AllowedIndicesBySpecies[newSpeciesIndex2_1][mirrorSymSigmaIndex1] && m_AllowedIndicesBySpecies[newSpeciesIndex2_2][mirrorSymSigmaIndex2];
    
      isTwoSwapAllowed = ((sigmaIndex1 != mirrorSymSigmaIndex1) && (sigmaIndex2 != mirrorSymSigmaIndex2) && (sigmaIndex1 != mirrorSymSigmaIndex2) && (sigmaIndex2 != mirrorSymSigmaIndex1) && (sigmaIndex1 != sigmaIndex2) && (mirrorSymSigmaIndex1 != mirrorSymSigmaIndex2));
      
      if((sigmaIndex1 <= m_MiddleStartLimit) || (sigmaIndex1 > m_MiddleEndLimit)){
          NotMiddleLayers = true;
      }
      else {
          NotMiddleLayers = false;
      }
      
      NotMiddleLayers &= ((sigmaIndex2 <= m_MiddleStartLimit)|| (sigmaIndex2 > m_MiddleEndLimit));
      NotMiddleLayers &= ((mirrorSymSigmaIndex1 <= m_MiddleStartLimit) || (mirrorSymSigmaIndex1 > m_MiddleEndLimit));
      NotMiddleLayers &= ((mirrorSymSigmaIndex2 <= m_MiddleStartLimit) || (mirrorSymSigmaIndex2 > m_MiddleEndLimit));
      
      isSwapAllowed &= NotMiddleLayers;
      isSwapAllowed &= isTwoSwapAllowed; 
    }        
    
    m_DoubleEvent.setSwap(sigmaIndex1, sigmaIndex2, mirrorSymSigmaIndex1, mirrorSymSigmaIndex2);
    return m_DoubleEvent;
  }

}

