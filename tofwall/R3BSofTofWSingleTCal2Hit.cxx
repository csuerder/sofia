// -----------------------------------------------------------------
// -----         R3BSofTofWSingleTCal2Hit source file          -----
// -----    Created 15/02/20  by J.L. Rodriguez-Sanchez        -----
// -----------------------------------------------------------------
#include "R3BSofTofWSingleTCal2Hit.h"
#include "R3BTGeoPar.h"

// R3BSofTofWSingleTCal2Hit: Default Constructor --------------------------
R3BSofTofWSingleTCal2Hit::R3BSofTofWSingleTCal2Hit()
    : FairTask("R3BSof TofW tcal2hit Task", 1)
    , fTCalDataCA(NULL)
    , fHitDataCA(NULL)
    , fExpId(467)
    , fOnline(kFALSE)
    , fTof_lise(43.)
{
}

// R3BSofTofWSingleTCal2Hit: Standard Constructor --------------------------
R3BSofTofWSingleTCal2Hit::R3BSofTofWSingleTCal2Hit(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fTCalDataCA(NULL)
    , fHitDataCA(NULL)
    , fExpId(467)
    , fOnline(kFALSE)
    , fTof_lise(43.)
{
}

// Virtual R3BSofTofWSingleTCal2Hit: Destructor
R3BSofTofWSingleTCal2Hit::~R3BSofTofWSingleTCal2Hit()
{
    LOG(INFO) << "R3BSofTofWSingleTCal2Hit: Delete instance";
    if (fTCalDataCA)
        delete fTCalDataCA;
    if (fHitDataCA)
        delete fHitDataCA;
}

void R3BSofTofWSingleTCal2Hit::SetParContainers()
{
    FairRuntimeDb* rtdb = FairRuntimeDb::instance();

    fTofWHitPar = (R3BSofTofWHitPar*)rtdb->getContainer("tofwHitPar");
    if (!fTofWHitPar)
    {
        LOG(ERROR) << "R3BSofTofWSingleTCal2Hit::SetParContainers() : Could not get access to tofwHitPar-Container.";
        return;
    }
    else
        LOG(INFO) << "R3BSofTofWSingleTCal2Hit::SetParContainers() : Container tofwHitPar found with "
                  << fTofWHitPar->GetNumSci() << " paddles.";

    fTofWGeoPar = (R3BTGeoPar*)rtdb->getContainer("tofwGeoPar");
    if (!fTofWGeoPar)
    {
        LOG(ERROR) << "R3BSofTofWSingleTCal2Hit::SetParContainers() : Could not get access to tofwGeoPar container.";
        return;
    }
    else
        LOG(INFO) << "R3BSofTofWSingleTCal2Hit::SetParContainers() : Container tofwGeoPar found.";
}

// -----   Public method Init   --------------------------------------------
InitStatus R3BSofTofWSingleTCal2Hit::Init()
{
    LOG(INFO) << "R3BSofTofWSingleTCal2Hit: Init";

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
    // fTofWHitPar->printParams();
    return kSUCCESS;
}

// -----   Public method ReInit   ----------------------------------------------
InitStatus R3BSofTofWSingleTCal2Hit::ReInit() { return kSUCCESS; }

// -----   Public method Execution   --------------------------------------------
void R3BSofTofWSingleTCal2Hit::Exec(Option_t* option)
{
    // At the moment we will use the expid to select the reconstruction
    // this should be changed in the future because expid is not necessary
    if (fExpId == 467 || fExpId == 444)
        S467();
    else if (fExpId == 455)
        S455();

    return;
}

// -----   Public method Experiment S455   --------------------------------------
void R3BSofTofWSingleTCal2Hit::S455()
{
    // Reset entries in output arrays, local arrays
    Reset();

    // Reading the Input -- SingleTCal Data
    Int_t nHits = fTCalDataCA->GetEntries();
    if (nHits < 1)
        return;

    // Data from cal level
    calDat = new R3BSofTofWSingleTcalData*[nHits];
    Int_t fPaddleId = 0; // from 1 to 28
    Double_t tofw = 0., posx = 0., posy = 0.;

    for (Int_t i = 0; i < nHits; i++)
    {
        calDat[i] = (R3BSofTofWSingleTcalData*)(fTCalDataCA->At(i));
        fPaddleId = calDat[i]->GetDetector();
        tofw = calDat[i]->GetRawTofNs();
        posy = calDat[i]->GetRawPosNs();

        posx = fTofWGeoPar->GetDimX() / 2.0 - 15. - (Double_t)(fPaddleId - 1) * 30.;
        posy = posy - fTofWHitPar->GetPosPar(fPaddleId);
        tofw = tofw - fTofWHitPar->GetTofPar(fPaddleId) + fTof_lise;

        AddHitData(fPaddleId, posx, posy, tofw);
    }
    if (calDat)
        delete calDat;
}

// -----   Public method Experiment S467   --------------------------------------
void R3BSofTofWSingleTCal2Hit::S467()
{
    // Reset entries in output arrays, local arrays
    Reset();

    // Reading the Input -- Cal Data --
    Int_t nHits = fTCalDataCA->GetEntries();
    if (nHits < 1)
        return;

    // Data from cal level
    calDat = new R3BSofTofWSingleTcalData*[nHits];
    Int_t fPaddleId = 0; // from 1 to 28
    Double_t tofw = 0., posx = 0., posy = 0.;
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
    {
        posx = fTofWGeoPar->GetDimX() / 2.0 - 15. - (Double_t)(fPaddleId - 1) * 30.; // x=0 at the gap of bars 14 and 15
        posy = posy - fTofWHitPar->GetPosPar(fPaddleId);
        tofw = tofw - fTofWHitPar->GetTofPar(fPaddleId) + fTof_lise;
        // TofPar is adjusted to align the time 0 with motor sweep runs,
        // And the Tof_lise is to adjust the difference of the flight path from sofsci to target setting-by-setting.
        AddHitData(fPaddleId, posx, posy, tofw);
    }
    if (calDat)
        delete calDat;
}

// -----   Public method Finish  ------------------------------------------------
void R3BSofTofWSingleTCal2Hit::Finish() {}

void R3BSofTofWSingleTCal2Hit::FinishEvent() {}

// -----   Public method Reset   ------------------------------------------------
void R3BSofTofWSingleTCal2Hit::Reset()
{
    LOG(DEBUG) << "Clearing TofWHitData structure";
    if (fHitDataCA)
        fHitDataCA->Clear();
}

// -----   Private method AddHitData  -------------------------------------------
R3BSofTofWHitData* R3BSofTofWSingleTCal2Hit::AddHitData(Int_t paddle, Double_t x, Double_t y, Double_t tof)
{
    // It fills the R3BSofTofWHitData
    TClonesArray& clref = *fHitDataCA;
    Int_t size = clref.GetEntriesFast();
    return new (clref[size]) R3BSofTofWHitData(paddle, x, y, tof);
}

ClassImp(R3BSofTofWSingleTCal2Hit)
