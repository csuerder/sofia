// ------------------------------------------------------------
// -----                  R3BSofFrsFillTree               -----
// -----              Fill SOFIA and FRS tree             ----
// -----           Created 07/08/20 by R. Taniuchi        -----
// ------------------------------------------------------------

/*
 * This task should fill histograms with SOFIA online data
 */

#include "R3BSofFrsFillTree.h"

R3BSofFrsFillTree::R3BSofFrsFillTree()
    : FairTask("SofFrsFillTree", 1)
    , fMappedItemsSci(NULL)
    , fTcalItemsSci(NULL)
    , fSingleTcalItemsSci(NULL)
    , fMusHitItems(NULL)
    , fMusCalItems(NULL)
    , fCalItemsMwpc0(NULL)
    , fNEvents(0)
    , fNbDetectors(4)
    , fNbChannels(3)
    , fIdS2(2)
    , fIdS8(3)
    , fIdCave(4)
{
}

R3BSofFrsFillTree::R3BSofFrsFillTree(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fMappedItemsSci(NULL)
    , fTcalItemsSci(NULL)
    , fSingleTcalItemsSci(NULL)
    , fMusHitItems(NULL)
    , fMusCalItems(NULL)
    , fCalItemsMwpc0(NULL)
    , fNEvents(0)
    , fNbDetectors(4)
    , fNbChannels(3)
    , fIdS2(2)
    , fIdS8(3)
    , fIdCave(4)
{
}
R3BSofFrsFillTree::~R3BSofFrsFillTree()
{
    LOG(INFO) << "R3BSofFrsFillTree::Delete instance";
    if (fMappedItemsSci)
        delete fMappedItemsSci;
    if (fTcalItemsSci)
        delete fTcalItemsSci;
    if (fSingleTcalItemsSci)
        delete fSingleTcalItemsSci;
    if (fMusHitItems)
        delete fMusHitItems;
    if (fTwimHitItems)
        delete fTwimHitItems;
    if (fMusCalItems)
        delete fMusCalItems;
    if (fCalItemsMwpc0)
        delete fCalItemsMwpc0;
}

InitStatus R3BSofFrsFillTree::Init()
{
    LOG(INFO) << "R3BSofFrsFillTree::Init ";

    // try to get a handle on the EventHeader. EventHeader may not be
    // present though and hence may be null. Take care when using.

    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(FATAL) << "R3BSofFrsFillTree::Init FairRootManager not found";
    header = (R3BEventHeader*)mgr->GetObject("R3BEventHeader");

    // Reading MusicCalPar from FairRuntimeDb
    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    if (!rtdb)
    {
        LOG(ERROR) << "FairRuntimeDb not opened!";
    }

    // --- ------------------------------------ --- //
    // --- get access to mapped data of the SCI --- //
    // --- ------------------------------------ --- //
    fFrsData = (TClonesArray*)mgr->GetObject("SofFrsData");
    if (!fFrsData)
    {
        return kFATAL;
    }

    // --- ------------------------------------ --- //
    // --- get access to mapped data of the SCI --- //
    // --- ------------------------------------ --- //
    fMappedItemsSci = (TClonesArray*)mgr->GetObject("SofSciMappedData");
    if (!fMappedItemsSci)
    {
        return kFATAL;
    }

    // --- ---------------------------------- --- //
    // --- get access to tcal data of the SCI --- //
    // --- ---------------------------------- --- //
    fTcalItemsSci = (TClonesArray*)mgr->GetObject("SofSciTcalData");
    if (!fTcalItemsSci)
    {
        return kFATAL;
    }

    // --- ----------------------------------------- --- //
    // --- get access to single tcal data of the SCI --- //
    // --- ----------------------------------------- --- //
    fSingleTcalItemsSci = (TClonesArray*)mgr->GetObject("SofSciSingleTcalData");
    if (!fSingleTcalItemsSci)
    {
        return kFATAL;
    }
    /*
    // get access to hit data of the MUSIC
    fMusHitItems = (TClonesArray*)mgr->GetObject("MusicHitData");
    if (!fMusHitItems)
        LOG(WARNING) << "R3BSofFrsFillTree: MusicHitData not found";

    // get access to cal data of the MUSIC
    fMusCalItems = (TClonesArray*)mgr->GetObject("MusicCalData");
    if (!fMusCalItems)
        LOG(WARNING) << "R3BSofFrsFillTree: MusicCalData not found";
    */

    // Twim
    fTwimHitItems = (TClonesArray*)mgr->GetObject("TwimHitData");
    if (!fTwimHitItems)
        LOG(WARNING) << "R3BSofFrsFillTree: TwimHitData not found";

    R3BSofTwimHitPar* fCal_TwimPar = (R3BSofTwimHitPar*)rtdb->getContainer("twimHitPar");
    if (!fCal_TwimPar)
    {
        LOG(ERROR) << "R3BSofTwimCal2HitPar::Init() Couldn't get handle on twimHitPar container";
    }
    //--- Parameter Container ---
    fNumSec = fCal_TwimPar->GetNumSec();        // Number of Sections
    fNumAnodes = fCal_TwimPar->GetNumAnodes();  // Number of anodes
    fNumParams = fCal_TwimPar->GetNumParZFit(); // Number of TwimParameters

    // Anodes that don't work set to zero
    TwimCalZParams = new TArrayF();
    Int_t array_size = fNumSec * fNumParams;
    TwimCalZParams->Set(array_size);
    TwimCalZParams = fCal_TwimPar->GetZHitPar(); // Array with the Cal parameters

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
    // Twim end

    // get access to cal data of the MWPC0
    fCalItemsMwpc0 = (TClonesArray*)mgr->GetObject("Mwpc0CalData");
    if (!fCalItemsMwpc0)
        LOG(WARNING) << "R3BSofFrsFillTree: Mwpc0CalData not found";

    // --- ------------------------------------ --- //
    // --- variables while looping over the data --- //
    // --- ------------------------------------ --- //
    // SofSci Mapped data
    multMapSci = new UChar_t[fNbDetectors * fNbChannels];
    // SofSci Tcal data
    iRawTimeNs = new Float_t[fNbDetectors * fNbChannels];

    // --- ------------------------------- --- //
    // ---    Create tree for detectors    --- //
    // --- ------------------------------- --- //
    // Add lines if necessary
    FrsTree = new TTree("FrsTree", "FrsTree");
    FrsTree->Branch("fNEvents", &fNEvents);
    FrsTree->Branch("tpat", &tpat);
    FrsTree->Branch("trigger", &trigger);
    FrsTree->Branch("MusicE", &MusicE);
    // FrsTree -> Branch("MusicDT", &MusicDT);
    FrsTree->Branch("MusicZ", &MusicZ);
    FrsTree->Branch("TwimE", &TwimE);
    FrsTree->Branch("TwimZ", &TwimZ);
    // FrsTree -> Branch("fNbDetectors", &fNbDetectors);
    // FrsTree -> Branch("fNbChannels", &fNbChannels);
    // FrsTree->Branch("multMapSci", multMapSci, "multMapSci[12]/b");
    // FrsTree->Branch("iRawTimeNs", iRawTimeNs, "iRawTimeNs[12]/F");
    FrsTree->Branch("xs2", &xs2);
    FrsTree->Branch("xpos", xpos, "xpos[3]/F");
    //
    FrsTree->Branch("TheBeta", &TheBeta);
    FrsTree->Branch("TheBrho", &TheBrho);
    FrsTree->Branch("TheAoQ", &TheAoQ);
    //
    FrsTree->Branch("Tof_wTref_S2_Cave", &Tof_wTref_S2_Cave);
    FrsTree->Branch("Beta_S2_Cave", &Beta_S2_Cave);
    FrsTree->Branch("Brho_S2_Cave", &Brho_S2_Cave);
    FrsTree->Branch("AoQ_S2_Cave", &AoQ_S2_Cave);
    FrsTree->Branch("MusicZ_S2_Cave", &MusicZ_S2_Cave);
    //
    FrsTree->Branch("Tof_wTref_S2_S8", &Tof_wTref_S2_S8);
    FrsTree->Branch("Beta_S2_S8", &Beta_S2_S8);
    FrsTree->Branch("Brho_S2_S8", &Brho_S2_S8);
    FrsTree->Branch("AoQ_S2_S8", &AoQ_S2_S8);
    FrsTree->Branch("MusicZ_S2_S8", &MusicZ_S2_S8);
    //
    FrsTree->Branch("Tof_wTref_S8_Cave", &Tof_wTref_S8_Cave);
    FrsTree->Branch("Beta_S8_Cave", &Beta_S8_Cave);
    FrsTree->Branch("Brho_S8_Cave", &Brho_S8_Cave);
    FrsTree->Branch("AoQ_S8_Cave", &AoQ_S8_Cave);
    FrsTree->Branch("MusicZ_S8_Cave", &MusicZ_S8_Cave);
    //
    return kSUCCESS;
}

void R3BSofFrsFillTree::Reset_Histo() { LOG(INFO) << "R3BSofFrsFillTree::Reset_Histo"; }

void R3BSofFrsFillTree::Exec(Option_t* option)
{
    fNEvents += 1;
    tpat = header->GetTpat();
    trigger = header->GetTrigger();
    if ((header->GetTpat() & 1) == 0 && header->GetTpat() != 0)
        return; // Phisics events should not be registered as tpat==0 but there are some..

    Int_t nHits;
    UShort_t iDet; // 0-bsed
    UShort_t iCh;  // 0-based

    // --- -------------- --- //
    // --- initialisation --- //
    // --- -------------- --- //
    for (UShort_t i = 0; i < fNbDetectors; i++)
    {
        for (UShort_t j = 0; j < fNbChannels; j++)
        {
            multMapSci[i * fNbChannels + j] = 0;
            iRawTimeNs[i * fNbChannels + j] = 0.;
        }
    }
    /*
    // --- -------------- --- //
    // --- MUSIC Hit data --- //
    // --- -------------- --- //
    if (fMusHitItems && fMusHitItems->GetEntriesFast() > 0)
      {
        nHits = fMusHitItems->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits; ihit++)
      {
            R3BMusicHitData* hit = (R3BMusicHitData*)fMusHitItems->At(ihit);
            if (!hit)
          continue;
            MusicE = hit->GetEave();
        //MusicZ = hit->GetZcharge();
      }
      }
    //
    */

    // --- -------------- --- //
    // --- Frs Ana data --- //
    // --- -------------- --- //
    if (fFrsData->GetEntriesFast() == 0)
        return;

    LOG(DEBUG) << "Entry" << fNEvents;

    nHits = fFrsData->GetEntriesFast();
    for (Int_t ihit = 0; ihit < nHits; ihit++)
    {
        R3BSofFrsData* hit = (R3BSofFrsData*)fFrsData->At(ihit);
        if (!hit)
            continue;
        // LOG(INFO) << "Sta: "<< hit->GetStaId() << " Sto: "<< hit->GetStoId() << " Z: "<<hit->GetZ()<<"
        // AoQ:"<<hit->GetAq();
        if (hit->GetStaId() == fIdS2)
        {
            if (hit->GetStoId() == fIdCave)
            {
                Beta_S2_Cave = hit->GetBeta();
                Brho_S2_Cave = hit->GetBrho();
                AoQ_S2_Cave = hit->GetAq();
                MusicZ_S2_Cave = hit->GetZ();
                if (xs2 < -1000.)
                    xs2 = hit->GetXS2();
            }
            else if (hit->GetStoId() == fIdS8)
            {
                Beta_S2_S8 = hit->GetBeta();
                Brho_S2_S8 = hit->GetBrho();
                AoQ_S2_S8 = hit->GetAq();
                MusicZ_S2_S8 = hit->GetZ();
                if (xs2 < -1000.)
                    xs2 = hit->GetXS2();
            }
            else
            {
                LOG(WARNING) << "Stop Sci not found. SciID=" << hit->GetStoId();
            }
        }
        else if (hit->GetStaId() == fIdS8 && hit->GetStoId() == fIdCave)
        {
            Beta_S8_Cave = hit->GetBeta();
            Brho_S8_Cave = hit->GetBrho();
            AoQ_S8_Cave = hit->GetAq();
            MusicZ_S8_Cave = hit->GetZ();
        }
        else
        {
            LOG(WARNING) << "Sci not found. StartSciID=" << hit->GetStaId() << ", StopSciID=" << hit->GetStoId();
        }
    } // End of the loop to obtain FrsData

    if (Beta_S2_Cave > 0.)
    {
        TheBeta = Beta_S2_Cave;
        TheBrho = Brho_S2_Cave;
        TheAoQ = AoQ_S2_Cave;
        MusicZ = MusicZ_S2_Cave;
    }
    else if (Beta_S2_S8 > 0.)
    {
        TheBeta = Beta_S2_S8;
        TheBrho = Brho_S2_S8;
        TheAoQ = AoQ_S2_S8;
        MusicZ = MusicZ_S2_S8;
    }
    else if (Beta_S8_Cave > 0.)
    {
        TheBeta = Beta_S8_Cave;
        TheBrho = Brho_S8_Cave;
        TheAoQ = AoQ_S8_Cave;
        MusicZ = MusicZ_S8_Cave;
    }

    if (fTwimHitItems && fTwimHitItems->GetEntriesFast() > 0)
    {
        nHits = fTwimHitItems->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BSofTwimHitData* Twimhit = (R3BSofTwimHitData*)fTwimHitItems->At(ihit);
            if (!Twimhit)
                continue;
            TwimE = Twimhit->GetEave();
        }
    }

    if (TheBeta > 0. && TwimE > 0.)
    {
        TwimZ = fTwimZ0 + fTwimZ1 * TMath::Sqrt(TwimE) * TheBeta + fTwimZ2 * TwimE * TheBeta * TheBeta;
    }

    // --- -------------- --- //
    // --- MUSIC Cal data --- //
    // --- -------------- --- //
    /*
    if (fMusCalItems && fMusCalItems->GetEntriesFast() > 0)
    {
        nHits = fMusCalItems->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BMusicCalData* hit = (R3BMusicCalData*)fMusCalItems->At(ihit);
            if (!hit)
                continue;
            if (hit->GetAnodeID() == 5)
                MusicDT = hit->GetDTime();
        }
    }*/

    if (fMappedItemsSci && fMappedItemsSci->GetEntriesFast() && fTcalItemsSci && fTcalItemsSci->GetEntriesFast())
    {
        // --- ------------------------- --- //
        // --- loop over sci mapped data --- //
        // --- ------------------------- --- //
        nHits = fMappedItemsSci->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BSofSciMappedData* hitmapped = (R3BSofSciMappedData*)fMappedItemsSci->At(ihit);
            if (!hitmapped)
                continue;
            iDet = hitmapped->GetDetector() - 1;
            iCh = hitmapped->GetPmt() - 1;
            multMapSci[iDet * fNbChannels + iCh]++;
            // fh1_finetime[iDet * fNbChannels + iCh]->Fill(hitmapped->GetTimeFine());
        }

        // --- ----------------------- --- //
        // --- loop over sci tcal data --- //
        // --- ----------------------- --- //
        if (fTcalItemsSci)
        {
            nHits = fTcalItemsSci->GetEntriesFast();
            for (Int_t ihit = 0; ihit < nHits; ihit++)
            {
                R3BSofSciTcalData* hittcal = (R3BSofSciTcalData*)fTcalItemsSci->At(ihit);
                if (!hittcal)
                    continue;
                iDet = hittcal->GetDetector() - 1;
                iCh = hittcal->GetPmt() - 1;
                iRawTimeNs[iDet * fNbChannels + iCh] = hittcal->GetRawTimeNs();
            }
        }
    }
    FrsTree->Fill();
}

// -----   Public method Reset   ------------------------------------------------
void R3BSofFrsFillTree::Reset() {}

// -----   Public method Finish   -----------------------------------------------
void R3BSofFrsFillTree::FinishEvent()
{
    for (UShort_t i = 0; i < fNbDetectors; i++)
    {
        for (UShort_t j = 0; j < fNbChannels; j++)
        {
            multMapSci[i * fNbChannels + j] = 0;
            iRawTimeNs[i * fNbChannels + j] = 0.;
        }
    }
    // Init branch values
    MusicZ = -5000., MusicE = -5000.;
    MusicDT = -5000.;
    TwimE = -5000., TwimZ = -5000.;
    xs2 = -5000.;
    xpos[0] = -5000.;
    xpos[1] = -5000.;
    xpos[2] = -5000.;
    MusicZ_S2_Cave = -5000.;
    Tof_wTref_S2_Cave = -5000., Beta_S2_Cave = -5000., Brho_S2_Cave = -5000.;
    MusicZ_S2_S8 = -5000.;
    Tof_wTref_S2_S8 = -5000., Beta_S2_S8 = -5000., Brho_S2_S8 = -5000.;
    MusicZ_S8_Cave = -5000.;
    Tof_wTref_S8_Cave = -5000., Beta_S8_Cave = -5000., Brho_S8_Cave = -5000.;
    AoQ_S2_Cave = -5000., AoQ_S2_S8 = -5000., AoQ_S8_Cave = -5000.;
    TheBeta = -5000., TheBrho = -5000., TheAoQ = -5000.;
    tpat = 0, trigger = 0;
    //
    if (fMappedItemsSci)
    {
        fMappedItemsSci->Clear();
    }
    if (fTcalItemsSci)
    {
        fTcalItemsSci->Clear();
    }
    if (fMusHitItems)
    {
        fMusHitItems->Clear();
    }
    if (fMusCalItems)
    {
        fMusCalItems->Clear();
    }
    if (fCalItemsMwpc0)
    {
        fCalItemsMwpc0->Clear();
    }
}

void R3BSofFrsFillTree::FinishTask() { FrsTree->Write(); }

ClassImp(R3BSofFrsFillTree)
