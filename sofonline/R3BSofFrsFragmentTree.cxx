// ------------------------------------------------------------
// -----                R3BSofFrsFragmentTree             -----
// -----              Fill FRS and SOFIA tree             -----
// -----           Created 07/08/20 by R. Taniuchi        -----
// ------------------------------------------------------------

/*
 * This task should fill histograms with SOFIA online data
 */

#include "R3BSofFrsFragmentTree.h"

R3BSofFrsFragmentTree::R3BSofFrsFragmentTree()
    : FairTask("SofFrsFragmentTree", 1)
    , fFrsData(NULL)
    , fFragData(NULL)
    //, fMappedItemsSci(NULL)
    //, fTcalItemsSci(NULL)
    //, fSingleTcalItemsSci(NULL)
    //, fMusHitItems(NULL)
    //, fMusCalItems(NULL)
    , fHitItemsMwpc0(NULL)
    , fHitItemsMwpc1(NULL)
    , fHitItemsMwpc2(NULL)
    , fHitItemsMwpc3(NULL)
    , fNEvents(0)
    //, fNbDetectors(4)
    //, fNbChannels(3)
    , fIdS2(2)
    , fIdS8(3)
    , fIdCave(4)
{
}

R3BSofFrsFragmentTree::R3BSofFrsFragmentTree(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fFrsData(NULL)
    , fFragData(NULL)
    //, fMappedItemsSci(NULL)
    //, fTcalItemsSci(NULL)
    //, fSingleTcalItemsSci(NULL)
    //, fMusHitItems(NULL)
    //, fMusCalItems(NULL)
    , fHitItemsMwpc0(NULL)
    , fHitItemsMwpc1(NULL)
    , fHitItemsMwpc2(NULL)
    , fHitItemsMwpc3(NULL)
    , fNEvents(0)
    //, fNbDetectors(4)
    //, fNbChannels(3)
    , fIdS2(2)
    , fIdS8(3)
    , fIdCave(4)
{
}
R3BSofFrsFragmentTree::~R3BSofFrsFragmentTree()
{
    LOG(INFO) << "R3BSofFrsFragmentTree::Delete instance";
    if (fFrsData)
        delete fFrsData;
    if (fFragData)
        delete fFragData;
    /*
    if (fMappedItemsSci)
        delete fMappedItemsSci;
    if (fTcalItemsSci)
        delete fTcalItemsSci;
    if (fSingleTcalItemsSci)
    delete fSingleTcalItemsSci;*/
    if (fMusHitItems)
      delete fMusHitItems;
    if (fTwimHitItems)
        delete fTwimHitItems;
    /*    if (fMusCalItems)
        delete fMusCalItems;
    */
    if (fTofWHitDataCA)
        delete fTofWHitDataCA;
    if (fHitItemsMwpc0)
        delete fHitItemsMwpc0;
    if (fHitItemsMwpc1)
        delete fHitItemsMwpc1;
    if (fHitItemsMwpc2)
        delete fHitItemsMwpc2;
    if (fHitItemsMwpc3)
        delete fHitItemsMwpc3;
}

InitStatus R3BSofFrsFragmentTree::Init()
{
    LOG(INFO) << "R3BSofFrsFragmentTree::Init ";

    // try to get a handle on the EventHeader. EventHeader may not be
    // present though and hence may be null. Take care when using.

    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(FATAL) << "R3BSofFrsFragmentTree::Init FairRootManager not found";
    header = (R3BEventHeader*)mgr->GetObject("R3BEventHeader");
    if (!header)
      header = (R3BEventHeader*)mgr->GetObject("EventHeader.");

    // Reading MusicCalPar from FairRuntimeDb
    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    if (!rtdb)
    {
        LOG(ERROR) << "FairRuntimeDb not opened!";
    }

    // --- ------------------------------------ --- //
    // --- get access to Ana data --- //
    // --- ------------------------------------ --- //
    fFrsData = (TClonesArray*)mgr->GetObject("SofFrsData");
    if (!fFrsData)
    {
        return kFATAL;
    }
    fFragData = (TClonesArray*)mgr->GetObject("SofTrackingData");
    if (!fFragData)
    {
        return kFATAL;
    }
    std::cout<<std::endl<<  __cplusplus << std::endl;
    /*
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
    */
    
    // get access to hit data of the MUSIC
    fMusHitItems = (TClonesArray*)mgr->GetObject("MusicHitData");
    if (!fMusHitItems)
        LOG(WARNING) << "R3BSofFrsFragmentTree: MusicHitData not found";
    /*
    // get access to cal data of the MUSIC
    fMusCalItems = (TClonesArray*)mgr->GetObject("MusicCalData");
    if (!fMusCalItems)
        LOG(WARNING) << "R3BSofFrsFragmentTree: MusicCalData not found";
    */

    // Twim
    fTwimHitItems = (TClonesArray*)mgr->GetObject("TwimHitData");
    if (!fTwimHitItems)
        LOG(WARNING) << "R3BSofFrsFragmentTree: TwimHitData not found";

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
    /*
    fCalItemsMwpc0 = (TClonesArray*)mgr->GetObject("Mwpc0CalData");
    if (!fCalItemsMwpc0)
        LOG(WARNING) << "R3BSofFrsFragmentTree: Mwpc0CalData not found";
    */
    fHitItemsMwpc0 = (TClonesArray*)mgr->GetObject("Mwpc0HitData");
    if (!fHitItemsMwpc0)
        LOG(WARNING) << "R3BSofFrsFragmentTree: Mwpc0HitData not found";

    fHitItemsMwpc1 = (TClonesArray*)mgr->GetObject("Mwpc1HitData");
    if (!fHitItemsMwpc1)
        LOG(WARNING) << "R3BSofFrsFragmentTree: Mwpc1HitData not found";

    fHitItemsMwpc2 = (TClonesArray*)mgr->GetObject("Mwpc2HitData");
    if (!fHitItemsMwpc2)
        LOG(WARNING) << "R3BSofFrsFragmentTree: Mwpc2HitData not found";

    fHitItemsMwpc3 = (TClonesArray*)mgr->GetObject("Mwpc3HitData");
    if (!fHitItemsMwpc3)
        LOG(WARNING) << "R3BSofFrsFragmentTree: Mwpc3HitData not found";

    fTofWHitDataCA = (TClonesArray*)mgr->GetObject("TofWHitData");
    if (!fTofWHitDataCA)
        LOG(WARNING) << "R3BSofFrsFragmentTree: TofWHitData not found";

    // --- ------------------------------------ --- //
    // --- variables while looping over the data --- //
    // --- ------------------------------------ --- //
    // SofSci Mapped data
    // multMapSci = new UChar_t[fNbDetectors * fNbChannels];
    // SofSci Tcal data
    // iRawTimeNs = new Float_t[fNbDetectors * fNbChannels];

    // --- ------------------------------- --- //
    // ---    Create tree for detectors    --- //
    // --- ------------------------------- --- //
    // Add lines if necessary
    Tree = new TTree("Tree", "Tree");
    Tree->Branch("fNEvents", &fNEvents);
    Tree->Branch("tpat", &tpat);
    Tree->Branch("trigger", &trigger);
    // Tree->Branch("MusicE", &MusicE);
    Tree->Branch("MusicZ", &MusicZ);
    Tree->Branch("MusicTheta", &MusicTheta);
    // Tree->Branch("TwimE", &TwimE);
    Tree->Branch("TwimZ", &TwimZ);
    Tree->Branch("TwimTheta", &TwimTheta);
    Tree->Branch("xs2", &xs2);
    Tree->Branch("xpos", xpos, "xpos[3]/F");
    //
    Tree->Branch("TheBeta", &TheBeta);
    Tree->Branch("TheGamma", &TheGamma);
    Tree->Branch("TheBrho", &TheBrho);
    Tree->Branch("TheAoQ", &TheAoQ);
    //
    Tree->Branch("Tof_wTref_S2_Cave", &Tof_wTref_S2_Cave);
    Tree->Branch("Beta_S2_Cave", &Beta_S2_Cave);
    Tree->Branch("Brho_S2_Cave", &Brho_S2_Cave);
    Tree->Branch("AoQ_S2_Cave", &AoQ_S2_Cave);
    Tree->Branch("MusicZ_S2_Cave", &MusicZ_S2_Cave);
    //
    Tree->Branch("Tof_wTref_S2_S8", &Tof_wTref_S2_S8);
    Tree->Branch("Beta_S2_S8", &Beta_S2_S8);
    Tree->Branch("Brho_S2_S8", &Brho_S2_S8);
    Tree->Branch("AoQ_S2_S8", &AoQ_S2_S8);
    Tree->Branch("MusicZ_S2_S8", &MusicZ_S2_S8);
    //
    Tree->Branch("Tof_wTref_S8_Cave", &Tof_wTref_S8_Cave);
    Tree->Branch("Beta_S8_Cave", &Beta_S8_Cave);
    Tree->Branch("Brho_S8_Cave", &Brho_S8_Cave);
    Tree->Branch("AoQ_S8_Cave", &AoQ_S8_Cave);
    Tree->Branch("MusicZ_S8_Cave", &MusicZ_S8_Cave);
    //
    Tree->Branch("Mw0_X", &Mw0_X);
    Tree->Branch("Mw0_Y", &Mw0_Y);
    Tree->Branch("Mw1_X", &Mw1_X);
    Tree->Branch("Mw1_Y", &Mw1_Y);
    Tree->Branch("Mw2_X", &Mw2_X);
    Tree->Branch("Mw2_Y", &Mw2_Y);
    Tree->Branch("Mw3_X", &Mw3_X);
    Tree->Branch("Mw3_Y", &Mw3_Y);
    // Tree->Branch("Tofw_Tof", &Tofw_Tof);
    Tree->Branch("Tofw_Paddle", &Tofw_Paddle);
    //
    Tree->Branch("FragZ", &FragZ);
    Tree->Branch("FragTof", &FragTof);
    Tree->Branch("FragBeta", &FragBeta);
    Tree->Branch("FragAoQ", &FragAoQ);
    Tree->Branch("FragBrho", &FragBrho);
    Tree->Branch("FragLength", &FragLength);
    return kSUCCESS;
}

void R3BSofFrsFragmentTree::Reset_Histo() { LOG(INFO) << "R3BSofFrsFragmentTree::Reset_Histo"; }

void R3BSofFrsFragmentTree::Exec(Option_t* option)
{
    fNEvents += 1;
    tpat = header->GetTpat();
    trigger = header->GetTrigger();
    if ((header->GetTpat() & 1) == 0 && header->GetTpat() != 0)
        return; // Phisics events should not be registered as tpat==0 but there are some..

    Int_t nHits;
    UShort_t iDet; // 0-bsed
    UShort_t iCh;  // 0-based
    /*
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
    */
    
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
	    MusicTheta = hit ->GetTheta();
        //MusicZ = hit->GetZcharge();
      }
      }
    //
    

    // --- -------------- --- //
    // --- Frs Ana data --- //
    // --- -------------- --- //
    if (fFrsData->GetEntriesFast() == 0)
    {
        return;
    }

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
    TheGamma = 1. / TMath::Sqrt(1 - TheBeta * TheBeta);

    if (fTwimHitItems && fTwimHitItems->GetEntriesFast() > 0)
    {
        nHits = fTwimHitItems->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BSofTwimHitData* Twimhit = (R3BSofTwimHitData*)fTwimHitItems->At(ihit);
            if (!Twimhit)
                continue;
            TwimE = Twimhit->GetEave();
	    TwimTheta = Twimhit ->GetTheta();
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
    }* /

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
    }*/

    if (fHitItemsMwpc0 && fHitItemsMwpc0->GetEntriesFast() > 0)
    {
        nHits = fHitItemsMwpc0->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            auto Mwpc0hit = (R3BSofMwpcHitData*)fHitItemsMwpc0->At(ihit);
            if (!Mwpc0hit)
                continue;
            Mw0_X = Mwpc0hit->GetX();
            Mw0_Y = Mwpc0hit->GetY();
        }
    }

    if (fHitItemsMwpc1 && fHitItemsMwpc1->GetEntriesFast() > 0)
    {
        nHits = fHitItemsMwpc1->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            auto Mwpc1hit = (R3BSofMwpcHitData*)fHitItemsMwpc1->At(ihit);
            if (!Mwpc1hit)
                continue;
            Mw1_X = Mwpc1hit->GetX();
            Mw1_Y = Mwpc1hit->GetY();
        }
    }

    if (fHitItemsMwpc2 && fHitItemsMwpc2->GetEntriesFast() > 0)
    {
        nHits = fHitItemsMwpc2->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            auto Mwpc2hit = (R3BSofMwpcHitData*)fHitItemsMwpc2->At(ihit);
            if (!Mwpc2hit)
                continue;
            Mw2_X = Mwpc2hit->GetX();
            Mw2_Y = Mwpc2hit->GetY();
        }
    }

    if (fHitItemsMwpc3 && fHitItemsMwpc3->GetEntriesFast() > 0)
    {
        nHits = fHitItemsMwpc3->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            auto Mwpc3hit = (R3BSofMwpcHitData*)fHitItemsMwpc3->At(ihit);
            if (!Mwpc3hit)
                continue;
            Mw3_X = Mwpc3hit->GetX();
            Mw3_Y = Mwpc3hit->GetY();
        }
    }
    /*
    if (fTofWHitDataCA && fTofWHitDataCA->GetEntriesFast() > 0)
    {
        nHits = fTofWHitDataCA->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            auto Tofwhit = (R3BSofTofWHitData*)fTofWHitDataCA->At(ihit);
            if (!Tofwhit)
                continue;
            FragTof = Tofwhit->GetTime();
            Tofw_Paddle = Tofwhit->GetPaddle();
        }
    }*/

    nHits = fFragData->GetEntriesFast();
    for (Int_t ihit = 0; ihit < nHits; ihit++)
    {
        R3BSofTrackingData* hit = (R3BSofTrackingData*)fFragData->At(ihit);
        if (!hit)
            continue;
        FragZ = hit->GetZ();
        FragAoQ = hit->GetAq();
        FragBeta = hit->GetBeta();
        FragBrho = hit->GetBrho();
        FragLength = hit->GetLength();
        FragTof = FragBrho / FragBeta;
        Tofw_Paddle = hit->GetPaddle();
    }
    // LOG(INFO) << MusicZ << " " << FragZ;
    Tree->Fill();
}

// -----   Public method Reset   ------------------------------------------------
void R3BSofFrsFragmentTree::Reset() {}

// -----   Public method Finish   -----------------------------------------------
void R3BSofFrsFragmentTree::FinishEvent()
{
    /*
      for (UShort_t i = 0; i < fNbDetectors; i++)
      {
          for (UShort_t j = 0; j < fNbChannels; j++)
          {
              multMapSci[i * fNbChannels + j] = 0;
              iRawTimeNs[i * fNbChannels + j] = 0.;
          }
      }
    */
    // Init branch values
    tpat = 0, trigger = 0;
    MusicZ = -5000., MusicE = -5000., MusicTheta = -5000.;
    MusicDT = -5000.;
    TwimE = -5000., TwimZ = -5000., TwimTheta = -5000.;
    xs2 = -5000.;
    xpos[0] = -5000.;
    xpos[1] = -5000.;
    xpos[2] = -5000.;
    Mw0_X = -5000., Mw0_Y = -5000.;
    Mw1_X = -5000., Mw1_Y = -5000.;
    Mw2_X = -5000., Mw2_Y = -5000.;
    Mw3_X = -5000., Mw3_Y = -5000.;
    MusicZ_S2_Cave = -5000.;
    Tof_wTref_S2_Cave = -5000., Beta_S2_Cave = -5000., Brho_S2_Cave = -5000.;
    MusicZ_S2_S8 = -5000.;
    Tof_wTref_S2_S8 = -5000., Beta_S2_S8 = -5000., Brho_S2_S8 = -5000.;
    MusicZ_S8_Cave = -5000.;
    Tof_wTref_S8_Cave = -5000., Beta_S8_Cave = -5000., Brho_S8_Cave = -5000.;
    AoQ_S2_Cave = -5000., AoQ_S2_S8 = -5000., AoQ_S8_Cave = -5000.;
    TheBeta = -5000., TheGamma = -5000., TheBrho = -5000., TheAoQ = -5000.;
    Tofw_Paddle = 0;
    FragZ = -5000., FragAoQ = -5000., FragTof = -5000., FragBeta = -5000., FragBrho = -5000., FragLength = -5000.;
    //
    /*
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
    }*/
    if (fHitItemsMwpc0)
    {
        fHitItemsMwpc0->Clear();
    }
    if (fHitItemsMwpc1)
    {
        fHitItemsMwpc1->Clear();
    }
    if (fHitItemsMwpc2)
    {
        fHitItemsMwpc2->Clear();
    }
    if (fHitItemsMwpc3)
    {
        fHitItemsMwpc3->Clear();
    }
}

void R3BSofFrsFragmentTree::FinishTask() { Tree->Write(); }

ClassImp(R3BSofFrsFragmentTree)
