// ----------------------------------------------------------------------
// -----                                                            -----
// -----                     R3BSofFragmentAnalysis                 -----
// -----             Created 09/02/20  by J.L. Rodriguez-Sanchez    -----
// ----------------------------------------------------------------------

#include "R3BSofFragmentAnalysis.h"

Double_t const c = 29.9792458;
TVector3 v1;

// R3BSofFragmentAnalysis: Default Constructor --------------------------
R3BSofFragmentAnalysis::R3BSofFragmentAnalysis()
    : FairTask("R3BSof Tracking Analysis", 1)
    //, fOffsetAq(0)
    //, fOffsetZ(0)
    //, fDist_mw3_tof(0.)
    //, fDist_start_glad(0.)
    , frho_Cave(0.)
    , fBfield_Glad(0.)
    , fTimeOffset(0.)
    , fTofWPos(560.) // 525: Centre
    , fMwpc0HitDataCA(NULL)
    , fMwpc1HitDataCA(NULL)
    , fMwpc2HitDataCA(NULL)
    , fMwpc3HitDataCA(NULL)
    , fTofWHitDataCA(NULL)
    , fTwimHitDataCA(NULL)
    , fTrackingDataCA(NULL)
    , fOnline(kFALSE)
{
}

// R3BSofFragmentAnalysisPar: Standard Constructor --------------------------
R3BSofFragmentAnalysis::R3BSofFragmentAnalysis(const TString& name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    //, fOffsetAq(0)
    //, fOffsetZ(0)
    //, fDist_mw3_tof(0.)
    //, fDist_start_glad(0.)
    , frho_Cave(0.)
    , fBfield_Glad(0.)
    , fTimeOffset(0.)
    , fTofWPos(560.) // 525: Centre, 560: 50Ca setting in s467
    , fMwpc0HitDataCA(NULL)
    , fMwpc1HitDataCA(NULL)
    , fMwpc2HitDataCA(NULL)
    , fMwpc3HitDataCA(NULL)
    , fTofWHitDataCA(NULL)
    , fTwimHitDataCA(NULL)
    , fTrackingDataCA(NULL)
    , fOnline(kFALSE)
{
}

// Virtual R3BSofFragmentAnalysis: Destructor
R3BSofFragmentAnalysis::~R3BSofFragmentAnalysis()
{
    LOG(INFO) << "R3BSofFragmentAnalysis: Delete instance";
    if (fMwpc0HitDataCA)
    {
        delete fMwpc0HitDataCA;
    }
    if (fMwpc1HitDataCA)
    {
        delete fMwpc1HitDataCA;
    }
    if (fMwpc2HitDataCA)
    {
        delete fMwpc2HitDataCA;
    }
    if (fMwpc3HitDataCA)
    {
        delete fMwpc3HitDataCA;
    }
    if (fTofWHitDataCA)
    {
        delete fTofWHitDataCA;
    }
    if (fTwimHitDataCA)
    {
        delete fTwimHitDataCA;
    }
    if (fTrackingDataCA)
    {
        delete fTrackingDataCA;
    }
}

void R3BSofFragmentAnalysis::SetParContainers()
{
    // Parameter Container
    // Reading softrackingAnaPar from FairRuntimeDb
    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    if (!rtdb)
    {
        LOG(ERROR) << "FairRuntimeDb not opened!";
    }
    //
    // Getting Twim Parameters
    R3BSofTwimHitPar* fTwimPar = (R3BSofTwimHitPar*)rtdb->getContainer("twimHitPar");
    if (!fTwimPar)
    {
        LOG(ERROR) << "R3BSofTwimCal2HitPar::Init() Couldn't get handle on twimHitPar container";
    }
    //--- Parameter Container ---
    Int_t fNumSec = fTwimPar->GetNumSec();        // Number of Sections
    Int_t fNumAnodes = fTwimPar->GetNumAnodes();  // Number of anodes
    Int_t fNumParams = fTwimPar->GetNumParZFit(); // Number of TwimParameters

    // Anodes that don't work set to zero
    TArrayF* TwimCalZParams = new TArrayF();
    Int_t array_size = fNumSec * fNumParams;
    TwimCalZParams->Set(array_size);
    TwimCalZParams = fTwimPar->GetZHitPar(); // Array with the Cal parameters

    // Parameters detector
    for (Int_t s = 0; s < fNumSec; s++)
        // Parameters detector
        if (fNumParams == 2)
        {
            fTwimZ0 = TwimCalZParams->GetAt(0);
            fTwimZ1 = TwimCalZParams->GetAt(1);
        }
        else if (fNumParams == 3)
        {
            fTwimZ0 = TwimCalZParams->GetAt(0);
            fTwimZ1 = TwimCalZParams->GetAt(1);
            fTwimZ2 = TwimCalZParams->GetAt(2);
        }
        else
            LOG(INFO) << "R3BSofTwimCal2Hit parameters for charge-Z cannot be used here, number of parameters: "
                      << fNumParams;
}

void R3BSofFragmentAnalysis::SetParameter()
{
    //--- Parameter Container ---
    // frho_Cave = 7.0;
    // fBfield_Glad = 4.0;
    fDist_mw3_tof = 72.0;
    fDist_start_glad = 65.5 + 163.4 + 118.;
    LOG(INFO) << "R3BSofFragmentAnalysis: Rho (Cave): " << frho_Cave;
    LOG(INFO) << "R3BSofFragmentAnalysis: B (Cave): " << fBfield_Glad;
}

// -----   Public method Init   --------------------------------------------
InitStatus R3BSofFragmentAnalysis::Init()
{
    LOG(INFO) << "R3BSofFragmentAnalysis: Init tracking analysis at Cave-C";

    // INPUT DATA
    FairRootManager* rootManager = FairRootManager::Instance();
    if (!rootManager)
    {
        return kFATAL;
    }

    fMwpc0HitDataCA = (TClonesArray*)rootManager->GetObject("Mwpc0HitData");
    if (!fMwpc0HitDataCA)
    {
        return kFATAL;
    }

    fMwpc1HitDataCA = (TClonesArray*)rootManager->GetObject("Mwpc1HitData");
    if (!fMwpc1HitDataCA)
    {
        return kFATAL;
    }

    fMwpc2HitDataCA = (TClonesArray*)rootManager->GetObject("Mwpc2HitData");
    if (!fMwpc2HitDataCA)
    {
        return kFATAL;
    }

    fMwpc3HitDataCA = (TClonesArray*)rootManager->GetObject("Mwpc3HitData");
    if (!fMwpc3HitDataCA)
    {
        return kFATAL;
    }

    fTwimHitDataCA = (TClonesArray*)rootManager->GetObject("TwimHitData");
    if (!fTwimHitDataCA)
    {
        return kFATAL;
    }

    fTofWHitDataCA = (TClonesArray*)rootManager->GetObject("TofWHitData");
    if (!fTofWHitDataCA)
    {
        return kFATAL;
    }

    // OUTPUT DATA
    fTrackingDataCA = new TClonesArray("R3BSofTrackingData", 2);
    if (!fOnline)
    {
        rootManager->Register("SofTrackingData", "GLAD Tracking Analysis", fTrackingDataCA, kTRUE);
    }
    else
    {
        rootManager->Register("SofTrackingData", "GLAD Tracking Analysis", fTrackingDataCA, kFALSE);
    }
    ReInit();
    SetParameter();
    return kSUCCESS;
}

// -----   Public method ReInit   ----------------------------------------------
InitStatus R3BSofFragmentAnalysis::ReInit()
{
    SetParContainers();
    return kSUCCESS;
}

// -----   Public method Execution   --------------------------------------------
void R3BSofFragmentAnalysis::Exec(Option_t* option)
{
    // Reset entries in output arrays, local arrays

    Double_t fZ = 0., fE = 0., fAq = 0.;
    Double_t Beta = 0., Brho_Cave = 0., Length = 0.;

    Int_t nHitMwpc0 = fMwpc0HitDataCA->GetEntries();
    Int_t nHitMwpc1 = fMwpc1HitDataCA->GetEntries();
    Int_t nHitMwpc2 = fMwpc2HitDataCA->GetEntries();
    Int_t nHitMwpc3 = fMwpc3HitDataCA->GetEntries();
    Int_t nHitTofW = fTofWHitDataCA->GetEntries();
    Int_t nHitTwim = fTwimHitDataCA->GetEntries();
    if (/*nHitMwpc0 < 1 || nHitMwpc1 < 1 || nHitMwpc2 < 1 ||*/ nHitMwpc3 < 1 || nHitTofW < 1 || nHitTwim < 1)
        return;
    // LOG(INFO) << "R3BSofFragmentAnalysis: nTwim: "<< nHitTwim << ", nTofW: " << nHitTofW << ", nMwpc: " << nHitMwpc ;

    R3BSofTofWHitData** HitTofW = new R3BSofTofWHitData*[nHitTofW];
    R3BSofTwimHitData** HitTwim = new R3BSofTwimHitData*[nHitTwim];
    R3BSofMwpcHitData** HitMwpc0 = new R3BSofMwpcHitData*[nHitMwpc0];
    R3BSofMwpcHitData** HitMwpc1 = new R3BSofMwpcHitData*[nHitMwpc1];
    R3BSofMwpcHitData** HitMwpc2 = new R3BSofMwpcHitData*[nHitMwpc2];
    R3BSofMwpcHitData** HitMwpc3 = new R3BSofMwpcHitData*[nHitMwpc3];

    // Added initial model
    Double_t ToF_Cave = 0.;
    Double_t mw3_x = 0., mw3_z = 0.;
    // Position at Cave-C behind GLAD: MWPC3
    for (Int_t i = 0; i < nHitMwpc3; i++)
    {
        HitMwpc3[i] = (R3BSofMwpcHitData*)(fMwpc3HitDataCA->At(i));
        mw3_x = HitMwpc3[i]->GetX() / 10. * cos(18. * TMath::DegToRad()) - 215.4;                        // cm
        mw3_z = 662. + HitMwpc3[i]->GetX() / 10. * sin(18. * TMath::DegToRad()) - 163.4 - fDist_mw3_tof; // cm
    }
    // Time from TofW
    for (Int_t i = 0; i < nHitTofW; i++)
    {
        HitTofW[i] = (R3BSofTofWHitData*)(fTofWHitDataCA->At(i));
        ToF_Cave = HitTofW[i]->GetTime();
        // std::cout <<" init: "<< HitTofW[i]->GetPaddle() << " "<< ToF_Cave << std::endl;
    }

    // std::cout << "R3BSofFragmentAnalysis: " << mw3_z << " " << mw3_x << " " << ToF_Cave ;//<< std::endl;

    if (ToF_Cave > 0. && mw3_z > 0. && mw3_x < 0.)
    {

        v1.SetXYZ(-1. * mw3_x, 0., mw3_z);

        double rho = 88.5969 / (2. * sin(v1.Theta() / 2.) * cos(14. * TMath::DegToRad() - v1.Theta() / 2.));
        Brho_Cave = 4.0 * 0.8 * rho * 0.01;
        Length = fDist_start_glad + sqrt(mw3_x * mw3_x + mw3_z * mw3_z) + fDist_mw3_tof;
        double vel = Length / ToF_Cave;
        Beta = vel / c;
        double gamma = 1. / sqrt(1. - Beta * Beta);
        fAq = Brho_Cave / (3.10716 * Beta * gamma);

        // std::cout << ToF_Cave << " " << vel << " " << Length << " " << fAq << " " << Brho_Cave << " " << Beta
        //          << std::endl;

        // Z from twim-music ------------------------------------
        Double_t countz = 0;
        for (Int_t i = 0; i < nHitTwim; i++)
        {
            HitTwim[i] = (R3BSofTwimHitData*)(fTwimHitDataCA->At(i));
            if (HitTwim[i]->GetZcharge() > 1)
            {
                // fZ = fZ + HitTwim[i]->GetZcharge();
                fE = fE + HitTwim[i]->GetEave();
                countz++;
            }
        }
        if (countz > 0)
        {
            fE = fE / countz;
            fZ = fTwimZ0 + fTwimZ1 * TMath::Sqrt(fE) * Beta + fTwimZ2 * fE * Beta * Beta;
        }
        else
        {
            fZ = 0.;
        }

        // Fill the data
        if (true) // if (fZ > 1 && fAq > 1. && Brho_Cave > 0. && Beta > 0.)
            AddData(fZ + fOffsetZ, fAq + fOffsetAq, Beta, Length, Brho_Cave);
    }

    if (HitTofW)
        delete HitTofW;
    if (HitMwpc0)
        delete HitMwpc0;
    if (HitMwpc1)
        delete HitMwpc1;
    if (HitMwpc2)
        delete HitMwpc2;
    if (HitMwpc3)
        delete HitMwpc3;
    if (HitTwim)
        delete HitTwim;
    return;
}

// -----   Protected method Finish   --------------------------------------------
void R3BSofFragmentAnalysis::Finish() {}

// -----   Public method Reset   ------------------------------------------------
void R3BSofFragmentAnalysis::Reset()
{
    LOG(DEBUG) << "Clearing SofTrackingData Structure";
    if (fTrackingDataCA)
        fTrackingDataCA->Clear();
}

// -----   Private method AddData  --------------------------------------------
R3BSofTrackingData* R3BSofFragmentAnalysis::AddData(Double_t z,
                                                    Double_t aq,
                                                    Double_t beta,
                                                    Double_t length,
                                                    Double_t brho)
{
    // It fills the R3BSofTrackingData
    TClonesArray& clref = *fTrackingDataCA;
    Int_t size = clref.GetEntriesFast();
    return new (clref[size]) R3BSofTrackingData(z, aq, beta, length, brho);
}
