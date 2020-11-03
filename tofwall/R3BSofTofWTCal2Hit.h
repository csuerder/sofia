// -----------------------------------------------------------------
// -----         R3BSofTofWTCal2Hit source file                -----
// -----    Created 15/02/20  by J.L. Rodriguez-Sanchez        -----
// -----------------------------------------------------------------

#ifndef R3BSofTofWTCal2Hit_H
#define R3BSofTofWTCal2Hit_H

#include "FairTask.h"
#include "R3BSofTofWTcalData.h"
#include "R3BSofTofWHitData.h"
#include "R3BSofTofWHitPar.h"
#include "TH1F.h"
#include <TRandom.h>

class TClonesArray;

class R3BSofTofWTCal2Hit : public FairTask
{

  public:
    /** Default constructor **/
    R3BSofTofWTCal2Hit();

    /** Standard constructor **/
    R3BSofTofWTCal2Hit(const char* name, Int_t iVerbose = 1);

    /** Destructor **/
    virtual ~R3BSofTofWTCal2Hit();

    /** Virtual method Exec **/
    virtual void Exec(Option_t* option);

    /** Virtual method Reset **/
    virtual void Reset();

    // Fair specific
    virtual void SetParContainers();

    /** Virtual method Init **/
    virtual InitStatus Init();

    /** Virtual method ReInit **/
    virtual InitStatus ReInit();

    /** Virtual method Finish **/
    virtual void Finish();

    void SetOnline(Bool_t option) { fOnline = option; }
    void SetTofLISE(Double_t tof) {Tof_lise = tof;}
    void SetTofWPos(Double_t pos) {TofWPosition = pos;}
    Double_t GetTofLISE() {return Tof_lise;}
    Double_t GetTofWPos() {return TofWPosition;}
        
  private:
    Bool_t fOnline; // Don't store data for online

    TClonesArray* fTCalDataCA; /**< Array with Cal input data. >*/
    TClonesArray* fHitDataCA;  /**< Array with Hit output data. >*/
    Double_t Tof_lise;// = 43.;
    Double_t TofWPosition;// = 560.
    R3BSofTofWHitPar* fTofWHitPar;
    
    /** Private method AddHitData **/
    // Adds a SofTofWHitData to the HitCollection
    R3BSofTofWHitData* AddHitData(Int_t paddle, Double_t x, Double_t y, Double_t tof);

  public:
    // Class definition
    ClassDef(R3BSofTofWTCal2Hit, 1)
};

#endif
