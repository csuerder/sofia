#include "R3BSofTofWMapped2Tcal.h"
#include "R3BSofTofWMappedData.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRunOnline.h"
#include "FairRuntimeDb.h"

R3BSofTofWMapped2Tcal::R3BSofTofWMapped2Tcal()
    : FairTask("R3BSofTofWMapped2Tcal", 1)
    , fMapped(NULL)
    , fTcalPar(NULL)
    , fTcal(new TClonesArray("R3BSofTofWTcalData"))
    , fNumTcal(0)
    , fOnline(kFALSE)
    , fNevent(0)
{
}

R3BSofTofWMapped2Tcal::~R3BSofTofWMapped2Tcal()
{
    if (fTcal)
    {
        delete fTcal;
    }
}

InitStatus R3BSofTofWMapped2Tcal::Init()
{

    LOG(INFO) << "R3BSofTofWMapped2Tcal: Init";

    FairRootManager* rm = FairRootManager::Instance();
    if (!rm)
    {
        LOG(ERROR) << "R3BSofTofWMapped2Tcal::Init() Couldn't instance the FairRootManager";
        return kFATAL;
    }

    // --- ----------------- --- //
    // --- INPUT MAPPED DATA --- //
    // --- ----------------- --- //

    // scintillator at S2 and cave C
    fMapped = (TClonesArray*)rm->GetObject("SofTofWMappedData"); // see Instance->Register() in R3BSofTofWReader.cxx
    if (!fMapped)
    {
        LOG(ERROR) << "R3BSofTofWMapped2Tcal::Init() Couldn't get handle on SofTofWMappedData container";
        return kFATAL;
    }
    else
        LOG(INFO) << " R3BSofTofWMapped2Tcal::Init() SofTofWMappedData items found";

    // --- -------------------------- --- //
    // --- CHECK THE TCALPAR VALIDITY --- //
    // --- -------------------------- --- //
    if (fTcalPar->GetNumSignals() == 0)
    {
        LOG(ERROR) << " R3BSofTofWMapped2Tcal::Init There are no Tcal parameters for SofTofW";
        return kFATAL;
    }
    else
    {
        LOG(INFO) << "R3BSofTofWMapped2Tcal::Init() : fNumSignals=" << fTcalPar->GetNumSignals();
        LOG(INFO) << " R3BSofTofWMapped2Tcal::Init() : fNumDetectors=" << fTcalPar->GetNumDetectors();
        LOG(INFO) << "  R3BSofTofWMapped2Tcal::Init() : fNumChannels=" << fTcalPar->GetNumChannels();
    }

    // --- ---------------- --- //
    // --- OUTPUT TCAL DATA --- //
    // --- ---------------- --- //

    // Register output array in tree
    if (!fOnline)
    {
        rm->Register("SofTofWTcalData", "SofTofW", fTcal, kTRUE);
    }
    else
    {
        rm->Register("SofTofWTcalData", "SofTofW", fTcal, kFALSE);
    }

    LOG(INFO) << "R3BSofTofWMapped2Tcal: Init DONE !";

    return kSUCCESS;
}

void R3BSofTofWMapped2Tcal::SetParContainers()
{
    fTcalPar = (R3BSofTcalPar*)FairRuntimeDb::instance()->getContainer("SofTofWTcalPar");
    if (!fTcalPar)
    {
        LOG(ERROR) << "R3BSofTofWMapped2Tcal::SetParContainers() : Could not get access to SofTofWTcalPar-Container.";
        return;
    }
    else
        LOG(INFO) << "R3BSofTofWMapped2Tcal::SetParContainers() : SofTofWTcalPar-Container found with "
                  << fTcalPar->GetNumSignals() << " signals";
}

InitStatus R3BSofTofWMapped2Tcal::ReInit()
{
    SetParContainers();
    return kSUCCESS;
}

void R3BSofTofWMapped2Tcal::Exec(Option_t* option)
{
    UShort_t iDet;
    UShort_t iCh;
    UInt_t iTf;
    UInt_t iTc;
    Double_t tns;

    Int_t nHitsPerEvent_SofTofW = fMapped->GetEntries();
    for (int ihit = 0; ihit < nHitsPerEvent_SofTofW; ihit++)
    {
        R3BSofTofWMappedData* hit = (R3BSofTofWMappedData*)fMapped->At(ihit);
        if (!hit)
            continue;
        iDet = hit->GetDetector();
        iCh = hit->GetPmt();
        iTf = hit->GetTimeFine();
        iTc = hit->GetTimeCoarse();
        if ((iDet < 1) || (iDet > fTcalPar->GetNumDetectors()))
        {
            LOG(INFO) << "R3BSofTofWMapped2Tcal::Exec() : In SofTofWMappedData, iDet = " << iDet
                      << "is out of range, item skipped ";
            continue;
        }
        if ((iCh < 1) || (iCh > fTcalPar->GetNumChannels()))
        {
            LOG(INFO) << "R3BSofTofWMapped2Tcal::Exec() : In SofTofWMappedData, iCh = " << iCh
                      << "is out of range, item skipped ";
            continue;
        }
        tns = CalculateTimeNs(iDet, iCh, iTf, iTc);
        new ((*fTcal)[fNumTcal++]) R3BSofTofWTcalData(iDet, iCh, tns);
    }

    ++fNevent;
}

void R3BSofTofWMapped2Tcal::FinishEvent()
{
    fTcal->Clear();
    fNumTcal = 0;
}

void R3BSofTofWMapped2Tcal::FinishTask() {}

Double_t R3BSofTofWMapped2Tcal::CalculateTimeNs(UShort_t iDet, UShort_t iCh, UInt_t iTf, UInt_t iTc)
{
    UInt_t rank = iTf + fTcalPar->GetNumTcalParsPerSignal() * ((iDet - 1) * fTcalPar->GetNumChannels() + (iCh - 1));
    Double_t iPar = (Double_t)fTcalPar->GetSignalTcalParams(rank);
    Double_t r = (Double_t)rand.Rndm() - 0.5;
    Double_t iTf_ns;
    Double_t iTc_ns = (Double_t)iTc * 5.;

    if (r < 0)
    {
        Double_t iParPrev = fTcalPar->GetSignalTcalParams(rank - 1);
        iTf_ns = iPar + r * (iPar - iParPrev);
    }
    else
    {
        Double_t iParNext = fTcalPar->GetSignalTcalParams(rank + 1);
        iTf_ns = iPar + r * (iParNext - iPar);
    }

    return (iTc_ns - iTf_ns);
}

ClassImp(R3BSofTofWMapped2Tcal)
