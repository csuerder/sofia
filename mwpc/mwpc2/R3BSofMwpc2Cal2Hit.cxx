// -------------------------------------------------------------------------
// -----         R3BSofMwpc2Cal2Hit source file                        -----
// -----             Created 11/10/19  by J.L. Rodriguez-Sanchez       -----
// -------------------------------------------------------------------------

// ROOT headers
#include "TClonesArray.h"
#include "TMath.h"

// Fair headers
#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"

#include <iomanip>

// MWPC headers
#include "R3BSofMwpc2Cal2Hit.h"
#include "R3BSofMwpcCalData.h"
#include "R3BSofMwpcHitData.h"

// R3BSofMwpc2Cal2Hit: Default Constructor --------------------------
R3BSofMwpc2Cal2Hit::R3BSofMwpc2Cal2Hit()
    : FairTask("R3B Hit-MWPC2 Task", 1)
    , fMwpcCalDataCA1(NULL)
    , fMwpcCalDataCA2(NULL)
    , fMwpcHitDataCA(NULL)
    , fwx(3.125)   // in mm
    , fwy(5.000)   // in mm
    , fSize(200.0) // in mm
    , fOnline(kFALSE)
{
}

// R3BSofMwpc2Cal2Hit: Standard Constructor --------------------------
R3BSofMwpc2Cal2Hit::R3BSofMwpc2Cal2Hit(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fMwpcCalDataCA1(NULL)
    , fMwpcCalDataCA2(NULL)
    , fMwpcHitDataCA(NULL)
    , fwx(3.125)   // in mm
    , fwy(5.000)   // in mm
    , fSize(200.0) // in mm
    , fOnline(kFALSE)
{
}

// Virtual R3BSofMwpc2Cal2Hit: Destructor
R3BSofMwpc2Cal2Hit::~R3BSofMwpc2Cal2Hit()
{
    LOG(INFO) << "R3BSofMwpc2Cal2Hit: Delete instance";
    if (fMwpcCalDataCA1)
        delete fMwpcCalDataCA1;
    if (fMwpcCalDataCA2)
        delete fMwpcCalDataCA2;
    if (fMwpcHitDataCA)
        delete fMwpcHitDataCA;
}

// -----   Public method Init   --------------------------------------------
InitStatus R3BSofMwpc2Cal2Hit::Init()
{
    LOG(INFO) << "R3BSofMwpc2Cal2Hit: Init";

    // INPUT DATA
    FairRootManager* rootManager = FairRootManager::Instance();
    if (!rootManager)
    {
        return kFATAL;
    }

    fMwpcCalDataCA1 = (TClonesArray*)rootManager->GetObject("Mwpc1CalData");
    if (!fMwpcCalDataCA1)
    {
        return kFATAL;
    }

    fMwpcCalDataCA2 = (TClonesArray*)rootManager->GetObject("Mwpc2CalData");
    if (!fMwpcCalDataCA2)
    {
        return kFATAL;
    }

    // OUTPUT DATA
    // Hit data
    fMwpcHitDataCA = new TClonesArray("R3BSofMwpcHitData", 10);

    if (!fOnline)
    {
        rootManager->Register("Mwpc2HitData", "MWPC2 Hit", fMwpcHitDataCA, kTRUE);
    }
    else
    {
        rootManager->Register("Mwpc2HitData", "MWPC2 Hit", fMwpcHitDataCA, kFALSE);
    }

    return kSUCCESS;
}

// -----   Public method ReInit   ----------------------------------------------
InitStatus R3BSofMwpc2Cal2Hit::ReInit() { return kSUCCESS; }

// -----   Public method Execution   --------------------------------------------
void R3BSofMwpc2Cal2Hit::Exec(Option_t* option)
{
    // Reset entries in output arrays, local arrays
    Reset();

    // Reading the Input -- Cal Data --
    Int_t nHits = fMwpcCalDataCA2->GetEntries();
    if (!nHits)
        return;

    // Data from cal level
    R3BSofMwpcCalData** calData;
    calData = new R3BSofMwpcCalData*[nHits];
    Int_t planeId;
    Int_t padId;
    Int_t padmx = -1, padmxu = -1, padmxd = -1, padmy = -1;
    Double_t q = 0, qmxu = 0, qmxd = 0, qleft = 0, qright = 0;
    Double_t qmx = 0, qmy = 0, qdown = 0, qup = 0;
    Double_t x = -1000., y = -1000.;

    for (Int_t i = 0; i < Mw2PadsX; i++)
        fx[i] = 0;
    for (Int_t i = 0; i < Mw2PadsY; i++)
        fy[i] = 0;

    for (Int_t i = 0; i < nHits; i++)
    {
        calData[i] = (R3BSofMwpcCalData*)(fMwpcCalDataCA2->At(i));
        planeId = calData[i]->GetPlane();
        padId = calData[i]->GetPad(); // From 0 to 63 for X down and up
        q = calData[i]->GetQ();

        // FIXME: in November this should be OK!
        if (planeId == 1 || planeId == 2)
            fx[padId] += q; // Xup+Xdown

        if (q > qmx && (planeId == 1 || planeId == 2))
        {
            qmx = q;
            padmx = padId;
        }
    }
    // Reading the Input -- Cal Data for Y --
    Int_t nHits1 = fMwpcCalDataCA1->GetEntries();
    if (!nHits1)
        return;

    // Data from cal level
    R3BSofMwpcCalData** calData1;
    calData1 = new R3BSofMwpcCalData*[nHits1];

    for (Int_t i = 0; i < nHits1; i++)
    {
        calData1[i] = (R3BSofMwpcCalData*)(fMwpcCalDataCA1->At(i));
        planeId = calData1[i]->GetPlane();
        padId = calData1[i]->GetPad(); // From 0 to 63 for X down and up
        q = calData1[i]->GetQ();

        if (planeId == 3)
            fy[padId] = q;
        if (q > qmy && planeId == 3)
        {
            qmy = q;
            padmy = padId;
        }
    }
    // Add Hit data ----
    if (padmx > 1 && padmy > 1 && padmx + 1 < Mw2PadsX && padmy + 1 < Mw2PadsY && qmx > 0 && qmy > 0)
    {
        // FIXME: in November this should be OK!
        // Obtain position X ----
        qleft = (Double_t)fx[padmx - 1];
        qright = (Double_t)fx[padmx + 1];
        if (qleft > 0 && qright > 0)
            x = GetPositionX(qmx, padmx, qleft, qright);
        // Obtain position Y ----
        qdown = fy[padmy - 1];
        qup = fy[padmy + 1];
        if (qdown > 0 && qup > 0)
            y = GetPositionY(qmy, padmy, qdown, qup);

        AddHitData(x, y);
    }

    if (calData)
        delete calData;
    if (calData1)
        delete calData1;
    return;
}

// -----   Protected method to obtain the position X ----------------------------
Double_t R3BSofMwpc2Cal2Hit::GetPositionX(Double_t qmax, Int_t padmax, Double_t qleft, Double_t qright)
{
    Double_t a3 = TMath::Pi() * fwx / (TMath::ACosH(0.5 * (TMath::Sqrt(qmax / qleft) + TMath::Sqrt(qmax / qright))));
    // Double_t a2 = gRandom->Uniform(-fwx / 2,fwx / 2);
    Double_t a2 = (a3 / TMath::Pi()) * TMath::ATanH((TMath::Sqrt(qmax / qleft) - TMath::Sqrt(qmax / qright)) /
                                                    (2 * TMath::SinH(TMath::Pi() * fwx / a3)));

    return (-1. * padmax * fwx + (fSize / 2) - (fwx / 2) - a2); // Left is positive and right negative
}

// -----   Protected method to obtain the position Y ----------------------------
Double_t R3BSofMwpc2Cal2Hit::GetPositionY(Double_t qmax, Int_t padmax, Double_t qdown, Double_t qup)
{
    Double_t a3 = TMath::Pi() * fwy / (TMath::ACosH(0.5 * (TMath::Sqrt(qmax / qdown) + TMath::Sqrt(qmax / qup))));
    // Double_t a2 = gRandom->Uniform(-fwy / 2, fwy / 2);
    Double_t a2 = (a3 / TMath::Pi()) * TMath::ATanH((TMath::Sqrt(qmax / qdown) - TMath::Sqrt(qmax / qup)) /
                                                    (2 * TMath::SinH(TMath::Pi() * fwy / a3)));

    return (padmax * fwy - (fSize / 2) + (fwy / 2) + a2);
}

// -----   Public method Finish  ------------------------------------------------
void R3BSofMwpc2Cal2Hit::Finish() {}

// -----   Public method Reset   ------------------------------------------------
void R3BSofMwpc2Cal2Hit::Reset()
{
    LOG(DEBUG) << "Clearing Mwpc2HitData Structure";
    if (fMwpcHitDataCA)
        fMwpcHitDataCA->Clear();
}

// -----   Private method AddHitData  --------------------------------------------
R3BSofMwpcHitData* R3BSofMwpc2Cal2Hit::AddHitData(Double_t x, Double_t y)
{
    // It fills the R3BSofMwpcHitData
    TClonesArray& clref = *fMwpcHitDataCA;
    Int_t size = clref.GetEntriesFast();
    return new (clref[size]) R3BSofMwpcHitData(x, y);
}

ClassImp(R3BSofMwpc2Cal2Hit)
