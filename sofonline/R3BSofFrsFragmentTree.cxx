// ------------------------------------------------------------
// -----                  R3BSofFrsFragmentTree           -----
// -----           Fill SOFIA online histograms           -----
// ------------------------------------------------------------

/*
 * This task should fill histograms with SOFIA online data
 */

#include "R3BSofFrsFragmentTree.h"

R3BSofFrsFragmentTree::R3BSofFrsFragmentTree()
    : FairTask("SofFrsFragmentTree", 1)
    , fMappedItemsSci(NULL)
    , fTcalItemsSci(NULL)
    , fSingleTcalItemsSci(NULL)
    , fMusHitItems(NULL)
    , fMusCalItems(NULL)
    , fCalItemsMwpc0(NULL)
    , fTofWHitData(NULL)
    , fNEvents(0)
    , fNbDetectors(1)
    , fNbChannels(3)
    , fIdS2(0)
    , fIdS8(0)
    , fBrho0(7.1175) // For 40Ca setting in s467
    , fS2SciCoef0(0.)
    , fS2SciCoef1(-5.8)
{
}

R3BSofFrsFragmentTree::R3BSofFrsFragmentTree(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fMappedItemsSci(NULL)
    , fTcalItemsSci(NULL)
    , fSingleTcalItemsSci(NULL)
    , fMusHitItems(NULL)
    , fMusCalItems(NULL)
    , fCalItemsMwpc0(NULL)
    , fTofWHitData(NULL)
    , fNEvents(0)
    , fNbDetectors(1)
    , fNbChannels(3)
    , fIdS2(0)
    , fIdS8(0)
    , fBrho0(7.1175) // For 40Ca setting in s467
    , fS2SciCoef0(0.)
    , fS2SciCoef1(-5.8)
{
}

R3BSofFrsFragmentTree::R3BSofFrsFragmentTree(Double_t brho, Double_t S2Sci0, Double_t S2Sci1)
    : FairTask("SofFrsFragmentTree", 1)
    , fMappedItemsSci(NULL)
    , fTcalItemsSci(NULL)
    , fSingleTcalItemsSci(NULL)
    , fMusHitItems(NULL)
    , fMusCalItems(NULL)
    , fCalItemsMwpc0(NULL)
    , fTofWHitData(NULL)
    , fNEvents(0)
    , fNbDetectors(1)
    , fNbChannels(3)
    , fIdS2(0)
    , fIdS8(0)
    , fBrho0(brho)
    , fS2SciCoef0(S2Sci0)
    , fS2SciCoef1(S2Sci1)
{
}

R3BSofFrsFragmentTree::R3BSofFrsFragmentTree(const char* name,
                                             Int_t iVerbose,
                                             Double_t brho,
                                             Double_t S2Sci0,
                                             Double_t S2Sci1)
    : FairTask(name, iVerbose)
    , fMappedItemsSci(NULL)
    , fTcalItemsSci(NULL)
    , fSingleTcalItemsSci(NULL)
    , fMusHitItems(NULL)
    , fMusCalItems(NULL)
    , fCalItemsMwpc0(NULL)
    , fTofWHitData(NULL)
    , fNEvents(0)
    , fNbDetectors(1)
    , fNbChannels(3)
    , fIdS2(0)
    , fIdS8(0)
    , fBrho0(brho)
    , fS2SciCoef0(S2Sci0)
    , fS2SciCoef1(S2Sci1)
{
}

R3BSofFrsFragmentTree::~R3BSofFrsFragmentTree()
{
    LOG(INFO) << "R3BSofFrsFragmentTree::Delete instance";
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
    if (fTofWHitData)
        delete fTofWHitData;
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

    // FairRunOnline* run = FairRunOnline::Instance();
    // run->GetHttpServer()->Register("", this);

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

    // get access to cal data of the MWPC0
    fCalItemsMwpc0 = (TClonesArray*)mgr->GetObject("Mwpc0CalData");
    if (!fCalItemsMwpc0)
        LOG(WARNING) << "R3BSofFrsFragmentTree: Mwpc0CalData not found";

    // get access to cal data of the MWPC1
    fCalItemsMwpc1 = (TClonesArray*)mgr->GetObject("Mwpc1CalData");
    if (!fCalItemsMwpc1)
        LOG(WARNING) << "R3BSofFrsFragmentTree: Mwpc1CalData not found";

    // get access to cal data of the MWPC2
    fCalItemsMwpc2 = (TClonesArray*)mgr->GetObject("Mwpc2CalData");
    if (!fCalItemsMwpc2)
        LOG(WARNING) << "R3BSofFrsFragmentTree: Mwpc2CalData not found";

    // get access to cal data of the MWPC3
    fCalItemsMwpc3 = (TClonesArray*)mgr->GetObject("Mwpc3CalData");
    if (!fCalItemsMwpc3)
        LOG(WARNING) << "R3BSofFrsFragmentTree: Mwpc3CalData not found";

    // get access to hit data of the MWPC3
    fHitItemsMwpc3 = (TClonesArray*)mgr->GetObject("Mwpc3HitData");
    if (!fHitItemsMwpc3)
        LOG(WARNING) << "R3BSofFrsFragmentTree: Mwpc3HitData not found";

    for (Int_t i = 0; i < 4; i++)
        MwpcPos[i] = new TVector3(0, 0, MwpcZpos[i]);

    // --- --------------------------------- --- //
    // --- get access to data of the MUSICs  --- //
    // --- --------------------------------- --- //
    // get access to hit data of the MUSIC
    fMusHitItems = (TClonesArray*)mgr->GetObject("MusicHitData");
    if (!fMusHitItems)
        LOG(WARNING) << "R3BSofFrsFragmentTree: MusicHitData not found";

    // get access to cal data of the MUSIC
    fMusCalItems = (TClonesArray*)mgr->GetObject("MusicCalData");
    if (!fMusCalItems)
        LOG(WARNING) << "R3BSofFrsFragmentTree: MusicCalData not found";

    // get access to hit data of the TwinMUSIC
    fTwimHitItems = (TClonesArray*)mgr->GetObject("TwimHitData");
    if (!fTwimHitItems)
        LOG(WARNING) << "R3BSofFrsFragmentTree: TwimHitData not found";

    // --- ----------------------------------- --- //
    // --- get access to tcal data of the TofW --- //
    // --- ----------------------------------- --- //
    fTcalItemsTofW = (TClonesArray*)mgr->GetObject("SofTofWTcalData");
    if (!fTcalItemsTofW)
    {
        return kFATAL;
    }

    // --- ------------------------------------------ --- //
    // --- get access to single tcal data of the TOfW --- //
    // --- ------------------------------------------ --- //
    fSingleTcalItemsTofW = (TClonesArray*)mgr->GetObject("SofTofWSingleTcalData");
    if (!fSingleTcalItemsTofW)
    {
        return kFATAL;
    }

    /*
    fTofWHitData = (TClonesArray*)mgr->GetObject("R3BSofTofWHitData");
    if (!fTofWHitData)
        LOG(WARNING) << "R3BSofFrsFragmentTree: TofWHitData not found";
    */

    // --- ----------------------------------- --- //
    // --- get access to ana data of fragments --- //
    // --- ----------------------------------- --- //

    fFragmentData = (TClonesArray*)mgr->GetObject("SofTrackingData");
    if (!fFragmentData)
    {
        return kFATAL;
    }

    ////////

    // Reading MusicCalPar from FairRuntimeDb
    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    if (!rtdb)
    {
        LOG(ERROR) << "FairRuntimeDb not opened!";
    }

    R3BMusicHitPar* fCal_Par; /// **< Parameter container. >* //
    fCal_Par = (R3BMusicHitPar*)rtdb->getContainer("musicHitPar");
    if (!fCal_Par)
    {
        LOG(ERROR) << "R3BMusicCal2HitPar::Init() Couldn't get handle on musicHitPar container";
    }
    else
    {
        LOG(INFO) << "R3BMusicCal2HitPar:: musicHitPar container open";
    }
    //--- Parameter Container ---
    fNumAnodes = fCal_Par->GetNumAnodes();  // Number of anodes
    fNumParams = fCal_Par->GetNumParZFit(); // Number of Parameters
    LOG(INFO) << "R3BMusicCal2Hit: Nb parameters for charge-Z: " << fNumParams;
    CalZParams = new TArrayF();
    CalZParams->Set(fNumParams);
    CalZParams = fCal_Par->GetZHitPar(); // Array with the Cal parameters

    // Parameters detector
    if (fNumParams == 2)
    {
        fZ0 = CalZParams->GetAt(0);
        fZ1 = CalZParams->GetAt(1);
    }
    else if (fNumParams == 3)
    {
        fZ0 = CalZParams->GetAt(0);
        fZ1 = CalZParams->GetAt(1);
        fZ2 = CalZParams->GetAt(2);
    }
    else
        LOG(INFO) << "R3BMusicCal2Hit parameters for charge-Z cannot be used here, number of parameters: "
                  << fNumParams;
    // getting parameters for R3BMUSIC end

    // Twim par
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
    Tree = new TTree("Tree", "Tree");
    Tree->Branch("fNEvents", &fNEvents);
    Tree->Branch("tpat", &tpat);
    Tree->Branch("trigger", &trigger);
    Tree->Branch("MusicE", &MusicE);
    // Tree -> Branch("MusicDT", &MusicDT);
    Tree->Branch("MusicZ", &MusicZ);
    Tree->Branch("TwimE", &TwimE);
    Tree->Branch("TwimZ", &TwimZ);
    // Tree -> Branch("fNbDetectors", &fNbDetectors);
    // Tree -> Branch("fNbChannels", &fNbChannels);
    Tree->Branch("multMapSci", multMapSci, "multMapSci[12]/b");
    Tree->Branch("iRawTimeNs", iRawTimeNs, "iRawTimeNs[12]/F");
    Tree->Branch("xs2", &xs2);
    Tree->Branch("xpos", xpos, "xpos[3]/F");
    //
    Tree->Branch("TheBeta", &TheBeta);
    Tree->Branch("TheBrho", &TheBrho);
    Tree->Branch("TheAoQ", &TheAoQ);
    //
    Tree->Branch("Tof_wTref_S2_Cave", &Tof_wTref_S2_Cave);
    Tree->Branch("Beta_S2_Cave", &Beta_S2_Cave);
    Tree->Branch("Brho_S2_Cave", &Brho_S2_Cave);
    Tree->Branch("AoQ_S2_Cave", &AoQ_S2_Cave);
    //
    Tree->Branch("Tof_wTref_S2_S8", &Tof_wTref_S2_S8);
    Tree->Branch("Beta_S2_S8", &Beta_S2_S8);
    Tree->Branch("Brho_S2_S8", &Brho_S2_S8);
    Tree->Branch("AoQ_S2_S8", &AoQ_S2_S8);
    //
    Tree->Branch("Tof_wTref_S8_Cave", &Tof_wTref_S8_Cave);
    Tree->Branch("Beta_S8_Cave", &Beta_S8_Cave);
    Tree->Branch("Brho_S8_Cave", &Brho_S8_Cave);
    Tree->Branch("AoQ_S8_Cave", &AoQ_S8_Cave);
    //
    Tree->Branch("MwpcPos", MwpcPos);
    // Tree -> Branch("TofW","R3BSofTofWHitData" ,&fTofWHitData);
    // Tree -> Branch("hitST", (&hitST), 300000);
    // Tree -> Branch("tofwhittcal", &tofwhittcal);
    Tree->Branch("TofWId", &TofWId);
    Tree->Branch("TofWX", &TofWX);
    Tree->Branch("TofWY", &TofWY);
    Tree->Branch("TofWT", &TofWT);
    //
    // FragmentData
    Tree->Branch("FragZ", &FragZ);
    Tree->Branch("FragAoQ", &FragAoQ);
    Tree->Branch("FragBeta", &FragBeta);
    Tree->Branch("FragLength", &FragLength);
    Tree->Branch("FragBrho", &FragBrho);
    return kSUCCESS;
}

void R3BSofFrsFragmentTree::Reset_Histo() { LOG(INFO) << "R3BSofFrsFragmentTree::Reset_Histo"; }

void R3BSofFrsFragmentTree::Exec(Option_t* option)
{
    fNEvents += 1;
    //
    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(FATAL) << "R3BSofFrsFragmentTree::Exec FairRootManager not found";

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

    ///
    // if (fCalItemsMwpc0 && fCalItemsMwpc1 && fCalItemsMwpc2 && fCalItemsMwpc3){
    // if(fCalItemsMwpc3->GetEntriesFast()!=1) continue;
    /*
    if (fHitItemsMwpc3 && fHitItemsMwpc3->GetEntriesFast() > 0)
      {
      R3BSofMwpcHitData* hit = (R3BSofMwpcCalData*)fHitItemsMwpc3->At(0);
      hit->GetX();
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
            // MusicZ = hit->GetZcharge();
        }
    }
    //

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

    if (MusicE < 0 && TwimE < 0)
        return; // End this event to reduce output root file size

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

    // --- -------------- --- //
    // --- TofW Hit data --- //
    // --- -------------- --- //
    /*
    if(fTofWHitData && fTofWHitData->GetEntriesFast() > 0){
      nHits = fTofWHitData->GetEntriesFast();
      tofhit = (R3BSofTofWHitData*)fTofWHitData->At(0);
      if(tofhit) std::cout << "TofW:" << nHits << ", Paddle:" << tofhit->GetTime() << std::endl;
    }
    */
    if (fSingleTcalItemsTofW && fSingleTcalItemsTofW->GetEntriesFast())
    {
        // --- ------------------------- --- //
        // --- loop over singe tcal data --- //
        // --- ------------------------- --- //
        nHits = fSingleTcalItemsTofW->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            // R3BSofTofWSingleTcalData * hitST = (R3BSofTofWSingleTcalData*)fSingleTcalItemsTofW->At(ihit);
            hitST.emplace_back((R3BSofTofWSingleTcalData*)fSingleTcalItemsTofW->At(ihit));
            // fh1_RawPos_AtSingleTcal[hitST->GetDetector()-1]->Fill(hitST->GetRawPosNs());
            // fh1_RawTof_AtSingleTcal[hitST->GetDetector()-1]->Fill(hitST->GetRawTofNs());
            // std::cout << "TofW:" << ihit << ", hitST_RawTof:" << hitST->GetRawTofNs() << std::endl;
            auto itr = hitST.end();
            --itr;
            TofWT.emplace_back((*itr)->GetRawTofNs());
            TofWY.emplace_back((*itr)->GetRawPosNs());
            // TofWX.emplace_back((Double_t)((*itr) -> GetDetector()));
            TofWId.emplace_back((*itr)->GetDetector());
            // std::cout << "SingleTcal:" << hitST.size() << " RawTof:" << (*itr)->GetRawTofNs()  << std::endl;
        } // end of loop over the singletcal data
    }
    if (fTcalItemsTofW && fTcalItemsTofW->GetEntriesFast())
    {
        // --- ------------------- --- //
        // --- loop over tcal data --- //
        // --- ------------------- --- //
        nHits = fTcalItemsTofW->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            // R3BSofTofWTcalData* hittcal = (R3BSofTofWTcalData*)fTcalItemsTofW->At(ihit);
            tofwhittcal.emplace_back((R3BSofTofWTcalData*)fTcalItemsTofW->At(ihit));
            // if (!hittcal)
            // continue;
            // if (*(hittcal.end())->GetPmt() == 3)
            //  continue;
            // iDet = hittcal->GetDetector() - 1;
            // iCh = hittcal->GetPmt() - 1;
            // iRawTimeNs[iDet * 2 + iCh] = hittcal->GetRawTimeNs();
            // std::cout << "TofW:" << ihit << ", hittcal_RawTime:" << hittcal->GetRawTimeNs() << std::endl;
        }
    }

    ///////////////////
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

        // --- ------------------------------ --- //
        // --- loop over sci single tcal data --- //
        // --- ------------------------------ --- //
        Int_t d, t;
        if (fSingleTcalItemsSci)
        {
            nHits = fSingleTcalItemsSci->GetEntriesFast();
            for (Int_t ihit = 0; ihit < nHits; ihit++)
            {
                R3BSofSciSingleTcalData* hitsingletcal = (R3BSofSciSingleTcalData*)fSingleTcalItemsSci->At(ihit);
                if (!hitsingletcal)
                    continue;
                d = hitsingletcal->GetDetector() - 1;
                if (d == fIdS2 - 1)
                    xs2 = hitsingletcal->GetRawPosNs() * fS2SciCoef1 + fS2SciCoef0;
                if (d > 0)
                    xpos[d - 1] = hitsingletcal->GetRawPosNs();
                if (fIdS2 > 0 && fIdS8 > 0)
                {
                    if (d == 2)
                        Tof_wTref_S2_S8 = hitsingletcal->GetRawTofNs_FromS2();
                    if (d == 3)
                    {
                        Tof_wTref_S2_Cave = hitsingletcal->GetRawTofNs_FromS2(); // same as toff
                        Tof_wTref_S8_Cave = hitsingletcal->GetRawTofNs_FromS8();
                    }
                }
            } // end of loop over the SingleTcalItems

            if (MusicE > 0)
            { // && xs2!=-5000.){ //&& toff!=-5000.) {
                if (Tof_wTref_S2_Cave > 0.)
                {
                    Beta_S2_Cave = 462.837731 / (Tof_wTref_S2_Cave - 1318.258541); // ToFCalib
                    Gamma_S2_Cave = 1. / (TMath::Sqrt(1. - (Beta_S2_Cave) * (Beta_S2_Cave)));
                    Brho_S2_Cave = fBrho0 * (1 + xs2 / 726.); //+mwpc0x/10./2000);
                    AoQ_S2_Cave = Brho_S2_Cave / (3.10716 * Gamma_S2_Cave * Beta_S2_Cave);
                }
                //
                if (Tof_wTref_S2_S8 > 0.)
                {
                    Beta_S2_S8 = 279.088230 / (Tof_wTref_S2_S8 - 694.519095); // ToFCalib
                    Gamma_S2_S8 = 1. / (TMath::Sqrt(1. - (Beta_S2_S8) * (Beta_S2_S8)));
                    Brho_S2_S8 = fBrho0 * (1 + xs2 / 726.); //+mwpc0x/10./2000);
                    AoQ_S2_S8 = Brho_S2_S8 / (3.10716 * Gamma_S2_S8 * Beta_S2_S8);
                }
                //
                if (Tof_wTref_S8_Cave > 0.)
                {
                    Beta_S8_Cave = 183.845298 / (Tof_wTref_S8_Cave - 623.62812); // ToFCalib
                    Gamma_S8_Cave = 1. / (TMath::Sqrt(1. - (Beta_S8_Cave) * (Beta_S8_Cave)));
                    Brho_S8_Cave = fBrho0 * (1 + (xs2 > -5000. ? xs2 : 0.) / 726.); //+mwpc0x/10./2000);
                    AoQ_S8_Cave = Brho_S8_Cave / (3.10716 * Gamma_S8_Cave * Beta_S8_Cave);
                }
                //
                if (Beta_S2_Cave > 0.)
                {
                    TheBeta = Beta_S2_Cave;
                    TheGamma = Gamma_S2_Cave;
                    TheBrho = Brho_S2_Cave;
                    TheAoQ = AoQ_S2_Cave;
                }
                else if (Beta_S2_S8 > 0.)
                {
                    TheBeta = Beta_S2_S8;
                    TheGamma = Gamma_S2_S8;
                    TheBrho = Brho_S2_S8;
                    TheAoQ = AoQ_S2_S8;
                }
                else if (Beta_S8_Cave > 0.)
                {
                    TheBeta = Beta_S8_Cave;
                    TheGamma = Gamma_S8_Cave;
                    TheBrho = Brho_S8_Cave;
                    TheAoQ = AoQ_S8_Cave;
                }
                if (TheBeta > 0.)
                {
                    MusicZ = fZ0 + fZ1 * TMath::Sqrt(MusicE) * TheBeta + fZ2 * MusicE * TheBeta * TheBeta;
                    TwimZ = fTwimZ0 + fTwimZ1 * TMath::Sqrt(TwimE) * TheBeta + fTwimZ2 * TwimE * TheBeta * TheBeta;
                }
                //
            } // if MusicE>0
        }
    }
    //
    if (fFragmentData && fFragmentData->GetEntriesFast() > 0)
    {
        nHits = fFragmentData->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BSofTrackingData* hit = (R3BSofTrackingData*)fFragmentData->At(ihit);
            if (!hit)
                continue;
            FragZ = hit->GetZ();
            FragAoQ = hit->GetAq();
            FragBeta = hit->GetBeta();
            FragLength = hit->GetLength();
            FragBrho = hit->GetBrho();
        }
    }
    //
    Tree->Fill();
}

// -----   Public method Reset   ------------------------------------------------
void R3BSofFrsFragmentTree::Reset() {}

// -----   Public method Finish   -----------------------------------------------
void R3BSofFrsFragmentTree::FinishEvent()
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
    // toff = -5000.;
    Tof_wTref_S2_Cave = -5000., Beta_S2_Cave = -5000., Gamma_S2_Cave = -5000., Brho_S2_Cave = -5000.;
    Tof_wTref_S2_S8 = -5000., Beta_S2_S8 = -5000., Gamma_S2_S8 = -5000., Brho_S2_S8 = -5000.;
    Tof_wTref_S8_Cave = -5000., Beta_S8_Cave = -5000., Gamma_S8_Cave = -5000., Brho_S8_Cave = -5000.;
    AoQ_S2_Cave = -5000., AoQ_S2_S8 = -5000., AoQ_S8_Cave = -5000.;
    TheBeta = -5000., TheGamma = -5000., TheBrho = -5000., TheAoQ = -5000.;
    tpat = 0, trigger = 0;
    //
    hitST.clear();
    tofwhittcal.clear();
    TofWX.clear();
    TofWId.clear();
    TofWY.clear();
    TofWT.clear();
    //
    FragZ = -5000., FragAoQ = -5000., FragBeta = -5000., FragLength = -5000., FragBrho = -5000.;
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
    if (fCalItemsMwpc1)
    {
        fCalItemsMwpc1->Clear();
    }
    if (fCalItemsMwpc2)
    {
        fCalItemsMwpc2->Clear();
    }
    if (fCalItemsMwpc3)
    {
        fCalItemsMwpc3->Clear();
    }
}

void R3BSofFrsFragmentTree::FinishTask() { Tree->Write(); }

ClassImp(R3BSofFrsFragmentTree)
