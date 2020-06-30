// ------------------------------------------------------------
// -----                R3BSofFrsFillTree            -----
// -----              Fill SOFIA and FRS tree             -----
// ------------------------------------------------------------

#ifndef R3BSofFrsFillTree_H
#define R3BSofFrsFillTree_H

#include "FairTask.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TMath.h"
#include <array>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

class TClonesArray;
class R3BEventHeader;

/**
 * This taks reads SCI data and plots online histograms
 */
class R3BSofFrsFillTree : public FairTask
{

  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BSofFrsFillTree();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BSofFrsFillTree(const char* name, Int_t iVerbose = 1);
    R3BSofFrsFillTree(Double_t brho = 7.1175, Double_t S2Sci0 = 0., Double_t S2Sci1 = -5.8);
    R3BSofFrsFillTree(const char* name, Int_t iVerbose =1, Double_t brho = 7.1175, Double_t S2Sci0 = 0., Double_t S2Sci1 = -5.8);
    
    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BSofFrsFillTree();

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

    void SetNbDetectors(Int_t ndets) {fNbDetectors = ndets;}
    void SetNbChannels(Int_t nchs) {fNbChannels = nchs;}
    void SetIdS2(Int_t id) {fIdS2 = id;}
    void SetIdS8(Int_t id) {fIdS8 = id;}
    void SetBrho0(Double_t brho) {fBrho0 = brho;}
    void SetS2Coef(Double_t S2Sci0, Double_t S2Sci1) {fS2SciCoef0 = S2Sci0; fS2SciCoef1 = S2Sci1;}
    Int_t GetNbDetectors() {return fNbDetectors;}
    Int_t GetNbChannels() {return fNbChannels;}
    Int_t GetIdS2() {return fIdS2;}
    Int_t GetIdS8() {return fIdS8;}
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
    TClonesArray* fTofwHitData;
    
    Int_t fNbDetectors;
    Int_t fNbChannels;
    Int_t fIdS2;
    Int_t fIdS8;
    Double_t fBrho0;  //Brho setting in FRS S2-S8
    Double_t fS2SciCoef0, fS2SciCoef1; //slope_calib = -5.8; // only for the s467, for S2 SofSci 


    Int_t fNumSec;
    Int_t fNumAnodes;
    Int_t fNumParams;
    Float_t fZ0 = 0., fZ1 = 0. , fZ2 = 0.; // CalibPar for R3BMUSIC
    TArrayF* CalZParams;
    Float_t fTwimZ0 = 0., fTwimZ1 = 0., fTwimZ2 = 0.; // CalibPar for Twim
    TArrayF* TwimCalZParams;

    // check for trigger should be done globablly (somewhere else)
    R3BEventHeader* header; /**< Event header.      */
    Int_t fNEvents;         /**< Event counter.     */

    // Tree
    TTree* FrsTree;
    Double_t MusicZ = -10000., MusicE = -10000.; //, MusicZ_betacorr = 0.;
    Double_t MusicDT = -1000000.;
    Double_t TwimE = -10000., TwimZ =10000.;
    UShort_t* multMapSci;
    Double_t* iRawTimeNs;
    Double_t xs2 = -10000.;
    //Double_t toff = -10000.;
    double Tof_wTref_S2_Cave = -10000., Beta_S2_Cave = -10000., Gamma_S2_Cave = -10000., Brho_S2_Cave = -10000.;
    double Tof_wTref_S2_S8 = -10000., Beta_S2_S8 = -10000., Gamma_S2_S8 = -10000., Brho_S2_S8 = -10000.;
    double Tof_wTref_S8_Cave = -10000., Beta_S8_Cave = -10000., Gamma_S8_Cave = -10000., Brho_S8_Cave = -10000.;
    double AoQ_S2_Cave = -10000., AoQ_S2_S8 = -10000., AoQ_S8_Cave = -10000.;
    
    /*
    // Canvas
    TCanvas** cSciMult;                // [fNbDetectors];
    TCanvas** cSciRawPos;              // [fNbDetectors];
    TCanvas** cMusicZvsRawPos;         // [fNbDetectors];
    TCanvas*  cMwpc0vsRawPos;
    TCanvas*  cMusicDTvsRawPos;
    TCanvas** cSciRawTof_FromS2;       // [fNbDetectors];
    TCanvas** cMusicZvsRawTof_FromS2;  // [fNbDetectors];
    TCanvas** cSciRawTof_FromS8;       // [fNbDetectors];
    TCanvas** cMusicZvsRawTof_FromS8;  // [fNbDetectors];
    TCanvas*  cMusicEvsBeta;
    TCanvas*  cTwimvsMusicZ_betacorrected;
    TCanvas*  cBeta_Correlation;
    TCanvas*  cAqvsx2;
    TCanvas*  cAqvsq;

    // Histograms for Mapped data : Fine Time and Mult
    TH1I** fh1_finetime;   // [fNbDetectors * NbChannels];
    TH2I** fh2_mult;       // [fNbDetectors];

    // Histograms for PosRaw Data at Tcal and SingleTcal
    TH1F** fh1_RawPos_AtTcalMult1;  // [fNbDetectors];
    TH1F** fh1_RawPos_AtSingleTcal; // [fNbDetectors];

    TH1D** fh1_RawTof_FromS2_AtTcalMult1;        // [fNbDetectors];
    TH1D** fh1_RawTof_FromS2_AtTcalMult1_wTref;  // [fNbDetectors];
    TH1D** fh1_RawTof_FromS2_AtSingleTcal_wTref; // [fNbDetectors];

    TH1D** fh1_RawTof_FromS8_AtTcalMult1;        // [fNbDetectors];
    TH1D** fh1_RawTof_FromS8_AtTcalMult1_wTref;  // [fNbDetectors];
    TH1D** fh1_RawTof_FromS8_AtSingleTcal_wTref; // [fNbDetectors];

    TH2F** fh2_Beta_Correlation;
    
    // Histogram for correlation with R3B-Music
    TH2F** fh2_MusZvsRawPos;          //[fNbDetectors];
    TH2F*  fh2_MusDTvsRawPos;
    TH2F** fh2_MusZvsRawTof_FromS2;   //[fNbDetectors];
    TH2F** fh2_MusZvsRawTof_FromS8;   //[fNbDetectors];
    TH2F*  fh2_TwimvsMusicZ_betacorrected;
    TH2F*  fh2_MusEvsBeta;
    TH2F*  fh2_Aqvsx2;
    TH2F*  fh2_Aqvsq;
    
    // Histogram for correlation with Mwpc0
    TH2F*  fh2_Mwpc0vsRawPos;
    */


  public:
    ClassDef(R3BSofFrsFillTree, 1)
};

#endif
