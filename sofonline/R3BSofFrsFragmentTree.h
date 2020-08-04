// ------------------------------------------------------------
// -----                R3BSofFrsFragmentTree            -----
// -----              Fill SOFIA and FRS tree             -----
// ------------------------------------------------------------

#ifndef R3BSofFrsFragmentTree_H
#define R3BSofFrsFragmentTree_H


#include "TCanvas.h"
#include "TClonesArray.h"
#include "TFolder.h"
#include "TH1.h"
#include "TH2.h"
#include "THttpServer.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TMath.h"
#include "TRandom.h"
#include "TTree.h"
#include "TVector3.h"

#include "FairTask.h"
#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRunOnline.h"
#include "FairRuntimeDb.h"

#include "R3BEventHeader.h"
#include "R3BMusicCalData.h"
#include "R3BMusicHitData.h"
#include "R3BMusicHitPar.h"
#include "R3BSofTwimCalData.h"
#include "R3BSofTwimHitData.h"
#include "R3BSofTwimHitPar.h"
#include "R3BSofMwpcCalData.h"
#include "R3BSofSciCalData.h"
#include "R3BSofSciMappedData.h"
#include "R3BSofSciSingleTcalData.h"
#include "R3BSofSciTcalData.h"
#include "R3BSofTofWHitData.h"
#include "R3BSofTofWSingleTcalData.h"
#include "R3BSofTofWTcalData.h"
#include "R3BSofFragmentAnalysis.h"

#include <array>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

class TClonesArray;
class R3BEventHeader;

/**
 * This taks reads SCI data and plots online histograms
 */
class R3BSofFrsFragmentTree : public FairTask
{

  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BSofFrsFragmentTree();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BSofFrsFragmentTree(const char* name, Int_t iVerbose = 1);
    R3BSofFrsFragmentTree(Double_t brho = 7.1175, Double_t S2Sci0 = 0., Double_t S2Sci1 = -5.8);
    R3BSofFrsFragmentTree(const char* name, Int_t iVerbose =1, Double_t brho = 7.1175, Double_t S2Sci0 = 0., Double_t S2Sci1 = -5.8);
    
    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BSofFrsFragmentTree();

    /**
     * Method for task initialization.
     * This function is called by the framework before
     * the event loop.
     * @return Initialization status. kSUCCESS, kERROR or kFATAL.
     */
    virtual InitStatus Init();

    /**
     * Method for event loop implementation.
     * Is called by the framework every time a new event is read.
     * @param option an execution option.
     */
    virtual void Exec(Option_t* option);

    /**
     * A method for finish of processing of an event.
     * Is called by the framework for each event after executing
     * the tasks.
     */
    virtual void FinishEvent();

    /**
     * Method for finish of the task execution.
     * Is called by the framework after processing the event loop.
     */
    virtual void FinishTask();

    /**
     * Methods to clean histograms.
     */
    virtual void Reset_Histo();

    /** Virtual method Reset **/
    virtual void Reset();

    void SetNbDetectors(UChar_t ndets) {fNbDetectors = ndets;}
    void SetNbChannels(UChar_t nchs) {fNbChannels = nchs;}
    void SetIdS2(UChar_t id) {fIdS2 = id;}
    void SetIdS8(UChar_t id) {fIdS8 = id;}
    void SetBrho0(Double_t brho) {fBrho0 = brho;}
    void SetS2Coef(Double_t S2Sci0, Double_t S2Sci1) {fS2SciCoef0 = S2Sci0; fS2SciCoef1 = S2Sci1;}
    UChar_t GetNbDetectors() {return fNbDetectors;}
    UChar_t GetNbChannels() {return fNbChannels;}
    UChar_t GetIdS2() {return fIdS2;}
    UChar_t GetIdS8() {return fIdS8;}
    Double_t GetBrho0() {return fBrho0;}
    Double_t GetS2Coef0() {return fS2SciCoef0;}
    Double_t GetS2Coef1() {return fS2SciCoef1;}

  private:
    TClonesArray* fMappedItemsSci;     /**< Array with mapped items. */
    TClonesArray* fTcalItemsSci;       /**< Array with tcal items. */
    TClonesArray* fSingleTcalItemsSci; /**< Array with tcal items. */
    TClonesArray* fMusHitItems;        /**< Array with MUSIC Hit items. */
    TClonesArray* fMusCalItems;        /**< Array with MUSIC Cal items. */
    TClonesArray* fTwimHitItems;        /**< Array with Twim Hit items. */
    TClonesArray* fCalItemsMwpc0;      /**< Array with cal items of mwpc0. */
    TClonesArray* fCalItemsMwpc1;      /**< Array with cal items of mwpc1. */
    TClonesArray* fCalItemsMwpc2;      /**< Array with cal items of mwpc2. */
    TClonesArray* fCalItemsMwpc3;      /**< Array with cal items of mwpc3. */
    TClonesArray* fHitItemsMwpc0;      /**< Array with cal items of mwpc0. */
    TClonesArray* fHitItemsMwpc1;      /**< Array with cal items of mwpc1. */
    TClonesArray* fHitItemsMwpc2;      /**< Array with cal items of mwpc2. */
    TClonesArray* fHitItemsMwpc3;      /**< Array with cal items of mwpc3. */
    TClonesArray* fTofWHitData;
    TClonesArray* fTcalItemsTofW;
    TClonesArray* fSingleTcalItemsTofW;
    TClonesArray* fFragmentData;
    
    UChar_t fNbDetectors;
    UChar_t fNbChannels;
    UChar_t fIdS2;
    UChar_t fIdS8;
    Double_t fBrho0;  //Brho setting in FRS S2-S8
    Double_t fS2SciCoef0, fS2SciCoef1; //slope_calib = -5.8; // only for the s467, for S2 SofSci 


    UChar_t fNumSec;
    UChar_t fNumAnodes;
    UChar_t fNumParams;
    Float_t fZ0 = 0., fZ1 = 0. , fZ2 = 0.; // CalibPar for R3BMUSIC
    TArrayF* CalZParams;
    Float_t fTwimZ0 = 0., fTwimZ1 = 0., fTwimZ2 = 0.; // CalibPar for Twim
    TArrayF* TwimCalZParams;

    // check for trigger should be done globablly (somewhere else)
    R3BEventHeader* header; /**< Event header.      */
    UInt_t fNEvents;         /**< Event counter.     */

    // Branch for FRS
    TTree* Tree;
    Float_t MusicZ = -5000., MusicE = -5000.; //, MusicZ_betacorr = 0.;
    Float_t MusicDT = -5000.;
    Float_t TwimE = -5000., TwimZ = -5000.;
    UChar_t* multMapSci;
    Float_t* iRawTimeNs;
    Float_t xs2 = -5000., xpos[3] = {-5000.};
    Float_t Tof_wTref_S2_Cave = -5000., Beta_S2_Cave = -5000., Gamma_S2_Cave = -5000., Brho_S2_Cave = -5000.;
    Float_t Tof_wTref_S2_S8 = -5000., Beta_S2_S8 = -5000., Gamma_S2_S8 = -5000., Brho_S2_S8 = -5000.;
    Float_t Tof_wTref_S8_Cave = -5000., Beta_S8_Cave = -5000., Gamma_S8_Cave = -5000., Brho_S8_Cave = -5000.;
    Float_t AoQ_S2_Cave = -5000., AoQ_S2_S8 = -5000., AoQ_S8_Cave = -5000.;
    Float_t TheBeta = -5000., TheGamma = -5000., TheBrho = -5000., TheAoQ = -5000.;
    UInt_t tpat = 0, trigger = 0;

    // Fragment branch
    // MWPC
    Double_t MwpcZpos[4] = {0.}; // Fix later
    TVector3 *MwpcPos[4];
    //TofW
    //R3BSofTofWHitData* tofhit;
    std::vector<R3BSofTofWSingleTcalData*> hitST;
    std::vector<R3BSofTofWTcalData*> tofwhittcal;
    std::vector<UChar_t> TofWId;
    std::vector<Float_t> TofWX;
    std::vector<Float_t> TofWY;
    std::vector<Float_t> TofWT;
    //
    //Frag ana
    Float_t FragZ = -5000., FragAoQ = -5000., FragBeta = -5000., FragLength = -5000., FragBrho = -5000.;


    
  public:
    ClassDef(R3BSofFrsFragmentTree, 1)
};

#endif
