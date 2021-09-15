// -----------------------------------------------------------------
// -----                                                       -----
// -----                R3BSofTwimMapped2CalPar                -----
// -----        Created 29/01/20  by J.L. Rodriguez-Sanchez    -----
// -----------------------------------------------------------------

#ifndef R3BSofTwimMapped2CalPar_H
#define R3BSofTwimMapped2CalPar_H

#include "FairTask.h"
#include "R3BSofTwimMapped2Cal.h"
#include "R3BSofTwimMappedData.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TH1F.h"

class TClonesArray;
class R3BSofTwimCalPar;

class R3BSofTwimMapped2CalPar : public FairTask
{
  public:
    /** Default constructor **/
    R3BSofTwimMapped2CalPar();

    /** Standard constructor **/
    R3BSofTwimMapped2CalPar(const TString& name,
                            Int_t iVerbose = 1,
                            const TString& namedeta = "Mwpc1",
                            const TString& namedetb = "Mwpc2");

    /** Destructor **/
    virtual ~R3BSofTwimMapped2CalPar();

    /** Virtual method Exec **/
    virtual void Exec(Option_t* option);

    /** Virtual method FinishEvent **/
    virtual void FinishEvent();

    /** Virtual method FinishTask **/
    virtual void FinishTask();

    /** Virtual method Reset **/
    virtual void Reset();

    /** Virtual method Init **/
    virtual InitStatus Init();

    /** Virtual method ReInit **/
    virtual InitStatus ReInit();

    /** Method to set up the position of MwpcA **/
    void SetPosMwpcA(Float_t pos) { fPosMwpcA = pos; }

    /** Method to set up the position of MwpcB **/
    void SetPosMwpcB(Float_t pos) { fPosMwpcB = pos; }

    /** Method to set up the position of Music **/
    void SetPosTwim(Float_t pos) { fPosTwim = pos; }

    /** Method to set up the limits for fit **/
    void SetFitLimits(Int_t left, Int_t right)
    {
        fLimit_left = left;
        fLimit_right = right;
	LOG(INFO)<<"R3BSofTwimMapped2CalPar Left"<<fLimit_left<<", Right"<<fLimit_right;
    }

  private:
    Int_t fNumSec;
    Int_t fNumAnodes;
    Int_t fMaxMult;
    Int_t fMinStatistics;
    Int_t fNumParams;
    Int_t fNumPosParams;
    Int_t fNumAnodesRef;
    Int_t fLimit_left;
    Int_t fLimit_right;
    Int_t fMaxSigma;
    TArrayF* CalParams;
    TArrayF* PosParams;

    TString fNameDetA;
    Float_t fPosMwpcA; // Position in the beam direction in mm
    TString fNameDetB;
    Float_t fPosMwpcB; // Position in the beam direction in mm
    Float_t fPosTwim;  // Position in the beam direction in mm

    Int_t mulanode[MAX_NB_TWIMSEC][MAX_NB_TWIMANODE + MAX_NB_TWIMTREF];
    Double_t fE[MAX_NB_TWIMSEC][MAX_MULT_TWIM_CAL][MAX_NB_TWIMANODE + MAX_NB_TWIMTREF];
    Double_t fDT[MAX_NB_TWIMSEC][MAX_MULT_TWIM_CAL][MAX_NB_TWIMANODE + MAX_NB_TWIMTREF];

    R3BSofTwimCalPar* fCal_Par;      /**< Parameter container. >*/
    TClonesArray* fTwimMappedDataCA; /**< Array with Twim Mapped-input data. >*/
    TClonesArray* fHitItemsMwpcA;    /**< Array with hit items. */
    TClonesArray* fHitItemsMwpcB;    /**< Array with hit items. */

    TGraph** fg_anode, **fg_anode_result;
    TGraph2D** fg_anode2d;

  public:
    // Class definition
    ClassDef(R3BSofTwimMapped2CalPar, 1)
};

#endif
