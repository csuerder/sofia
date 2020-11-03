// -----------------------------------------------------------------
// -----         R3BSofTofWTCal2Hit source file                -----
// -----    Created 15/02/20  by J.L. Rodriguez-Sanchez        -----
// -----------------------------------------------------------------

// ROOT headers
#include "TClonesArray.h"
#include "TMath.h"

// Fair headers
#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"

#include <iomanip>

// TofW headers
#include "R3BSofTofWHitData.h"
#include "R3BSofTofWSingleTcalData.h"
#include "R3BSofTofWTCal2Hit.h"

// R3BSofTofWTCal2Hit: Default Constructor --------------------------
R3BSofTofWTCal2Hit::R3BSofTofWTCal2Hit()
    : FairTask("R3BSof TofW tcal2hit Task", 1)
    , fTCalDataCA(NULL)
    , fHitDataCA(NULL)
    , fOnline(kFALSE)
    , Tof_lise(43.)
    , TofWPosition(560.)
{
}

// R3BSofTofWTCal2Hit: Standard Constructor --------------------------
R3BSofTofWTCal2Hit::R3BSofTofWTCal2Hit(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fTCalDataCA(NULL)
    , fHitDataCA(NULL)
    , fOnline(kFALSE)
    , Tof_lise(43.)
    , TofWPosition(560.)
{
}

// Virtual R3BSofTofWTCal2Hit: Destructor
R3BSofTofWTCal2Hit::~R3BSofTofWTCal2Hit()
{
    LOG(INFO) << "R3BSofTofWTCal2Hit: Delete instance";
    if (fTCalDataCA)
        delete fTCalDataCA;
    if (fHitDataCA)
        delete fHitDataCA;
}

void R3BSofTofWTCal2Hit::SetParContainers()
{
    fTofWHitPar = (R3BSofTofWHitPar*)FairRuntimeDb::instance()->getContainer("tofwHitPar");
    if (!fTofWHitPar)
    {
        LOG(ERROR) << "R3BSofTofWTCal2Hit::SetParContainers() : Could not get access to tofwHitPar-Container.";
        return;
    }
    else
        LOG(INFO) << "R3BSofTofWTCal2Hit::SetParContainers() : tofwHitPar-Container found with "
                  << fTofWHitPar->GetNumSci() << " paddles.";
    fTofWHitPar->printParams();
}

// -----   Public method Init   --------------------------------------------
InitStatus R3BSofTofWTCal2Hit::Init()
{
    LOG(INFO) << "R3BSofTofWTCal2Hit: Init";

    // INPUT DATA
    FairRootManager* rootManager = FairRootManager::Instance();
    if (!rootManager)
    {
        return kFATAL;
    }

    fTCalDataCA = (TClonesArray*)rootManager->GetObject("SofTofWSingleTcalData");
    if (!fTCalDataCA)
    {
        return kFATAL;
    }

    // OUTPUT DATA
    // Hit data
    fHitDataCA = new TClonesArray("R3BSofTofWHitData", 10);

    if (!fOnline)
    {
        rootManager->Register("TofWHitData", "TofW-Hit", fHitDataCA, kTRUE);
    }
    else
    {
        rootManager->Register("TofWHitData", "TofW-Hit", fHitDataCA, kFALSE);
    }

    return kSUCCESS;
}

// -----   Public method ReInit   ----------------------------------------------
InitStatus R3BSofTofWTCal2Hit::ReInit() { return kSUCCESS; }

// -----   Public method Execution   --------------------------------------------
void R3BSofTofWTCal2Hit::Exec(Option_t* option)
{
    // Reset entries in output arrays, local arrays
    Reset();

    // Reading the Input -- Cal Data --
    Int_t nHits = fTCalDataCA->GetEntries();
    if (nHits < 1)
        return;

    // Data from cal level
    R3BSofTofWSingleTcalData** calDat;
    calDat = new R3BSofTofWSingleTcalData*[nHits];
    Int_t fPaddleId = 0; // from 1 to 28
    Double_t tofw = 0., posy = 0.;
    Int_t mult = 0;

    for (Int_t i = 0; i < nHits; i++)
    {
        calDat[i] = (R3BSofTofWSingleTcalData*)(fTCalDataCA->At(i));
        fPaddleId = calDat[i]->GetDetector();
        if (fTofWHitPar->GetInUse(fPaddleId) != 1)
            continue;
        mult++;
        tofw = calDat[i]->GetRawTofNs();
        posy = calDat[i]->GetRawPosNs();
    }

    if (mult == 1)
        AddHitData(
            fPaddleId,
            (TofWPosition - 525.) - 14. * 30. + (Double_t)(fPaddleId - 1) * 30., // x=0 at the gap of bars 14 and 15
            posy - fTofWHitPar->GetPosPar(fPaddleId),
            tofw - fTofWHitPar->GetTofPar((Int_t)fPaddleId) + Tof_lise);

    if (calDat)
        delete calDat;
    return;
}

// -----   Public method Finish  ------------------------------------------------
void R3BSofTofWTCal2Hit::Finish() {}

// -----   Public method Reset   ------------------------------------------------
void R3BSofTofWTCal2Hit::Reset()
{
    LOG(DEBUG) << "Clearing TofWHitData structure";
    if (fHitDataCA)
        fHitDataCA->Clear();
}

// -----   Private method AddHitData  --------------------------------------------
R3BSofTofWHitData* R3BSofTofWTCal2Hit::AddHitData(Int_t paddle, Double_t x, Double_t y, Double_t tof)
{
    // It fills the R3BSofTofWHitData
    TClonesArray& clref = *fHitDataCA;
    Int_t size = clref.GetEntriesFast();
    return new (clref[size]) R3BSofTofWHitData(paddle, x, y, tof);
}

ClassImp(R3BSofTofWTCal2Hit)
