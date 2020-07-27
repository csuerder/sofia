// ------------------------------------------------------------
// -----                  R3BSofFrsFillTree           -----
// -----           Fill SOFIA online histograms           -----
// ------------------------------------------------------------

/*
 * This task should fill histograms with SOFIA online data
 */

#include "R3BSofFrsFillTree.h"
#include "R3BEventHeader.h"
#include "R3BMusicCalData.h"
#include "R3BMusicHitData.h"
#include "R3BMusicHitPar.h"
#include "R3BSofTwimCalData.h"
#include "R3BSofTwimHitData.h"
#include "R3BSofTwimHitPar.h"
#include "R3BSofMwpcCalData.h"
#include "R3BSofSciCalData.h"
#include "R3BSofSciMappedData.h"
#include "R3BSofSciSingleTcalData.h"
#include "R3BSofSciTcalData.h"
#include "THttpServer.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRunOnline.h"
#include "FairRuntimeDb.h"
#include "TCanvas.h"
#include "TFolder.h"
#include "TH1.h"
#include "TH2.h"
#include "TVector3.h"

#include "TClonesArray.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TMath.h"
#include "TRandom.h"
#include <array>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>

R3BSofFrsFillTree::R3BSofFrsFillTree()
    : FairTask("SofFrsFillTree", 1)
    , fMappedItemsSci(NULL)
    , fTcalItemsSci(NULL)
    , fSingleTcalItemsSci(NULL)
    , fMusHitItems(NULL)
    , fMusCalItems(NULL)
    , fCalItemsMwpc0(NULL)
    , fTofwHitData(NULL)
    , fNEvents(0)
    , fNbDetectors(1)
    , fNbChannels(3)
    , fIdS2(0)
    , fIdS8(0)
    , fBrho0(7.1175) //For 40Ca setting in s467
    , fS2SciCoef0(0.)
    , fS2SciCoef1(-5.8)
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
    , fTofwHitData(NULL)
    , fNEvents(0)
    , fNbDetectors(1)
    , fNbChannels(3)
    , fIdS2(0)
    , fIdS8(0)
    , fBrho0(7.1175) //For 40Ca setting in s467
    , fS2SciCoef0(0.)
    , fS2SciCoef1(-5.8)
{
}

R3BSofFrsFillTree::R3BSofFrsFillTree(Double_t brho, Double_t S2Sci0, Double_t S2Sci1)
    : FairTask("SofFrsFillTree", 1)
    , fMappedItemsSci(NULL)
    , fTcalItemsSci(NULL)
    , fSingleTcalItemsSci(NULL)
    , fMusHitItems(NULL)
    , fMusCalItems(NULL)
    , fCalItemsMwpc0(NULL)
    , fTofwHitData(NULL)
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

R3BSofFrsFillTree::R3BSofFrsFillTree(const char* name, Int_t iVerbose, Double_t brho, Double_t S2Sci0, Double_t S2Sci1)
    : FairTask(name, iVerbose)
    , fMappedItemsSci(NULL)
    , fTcalItemsSci(NULL)
    , fSingleTcalItemsSci(NULL)
    , fMusHitItems(NULL)
    , fMusCalItems(NULL)
    , fCalItemsMwpc0(NULL)
    , fTofwHitData(NULL)
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
    if (fTofwHitData)
        delete fTofwHitData;
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

    //FairRunOnline* run = FairRunOnline::Instance();
    //run->GetHttpServer()->Register("", this);

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

    // get access to hit data of the MUSIC
    fMusHitItems = (TClonesArray*)mgr->GetObject("MusicHitData");
    if (!fMusHitItems)
        LOG(WARNING) << "R3BSofFrsFillTree: MusicHitData not found";

    // get access to cal data of the MUSIC
    fMusCalItems = (TClonesArray*)mgr->GetObject("MusicCalData");
    if (!fMusCalItems)
        LOG(WARNING) << "R3BSofFrsFillTree: MusicCalData not found";

    // Reading MusicCalPar from FairRuntimeDb
    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    if (!rtdb)
    {
        LOG(ERROR) << "FairRuntimeDb not opened!";
    }

    R3BMusicHitPar* fCal_Par;      /// **< Parameter container. >* //
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
    fNumAnodes = fCal_Par->GetNumAnodes(); // Number of anodes
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
    FrsTree = new TTree("FrsTree", "FrsTree");
    FrsTree -> Branch("fNEvents", &fNEvents);
    FrsTree -> Branch("tpat", &tpat);
    FrsTree -> Branch("trigger", &trigger);
    FrsTree -> Branch("MusicE", &MusicE);
    //FrsTree -> Branch("MusicDT", &MusicDT);
    FrsTree -> Branch("MusicZ", &MusicZ);
    FrsTree -> Branch("TwimE", &TwimE);
    FrsTree -> Branch("TwimZ", &TwimZ);
    //FrsTree -> Branch("fNbDetectors", &fNbDetectors);
    //FrsTree -> Branch("fNbChannels", &fNbChannels);
    FrsTree -> Branch("multMapSci", multMapSci, "multMapSci[12]/b");
    FrsTree -> Branch("iRawTimeNs", iRawTimeNs, "iRawTimeNs[12]/F");
    FrsTree -> Branch("xs2", &xs2);
    FrsTree -> Branch("xpos", xpos, "xpos[3]/F");
    //
    FrsTree -> Branch("TheBeta", &TheBeta);
    FrsTree -> Branch("TheBrho", &TheBrho);
    FrsTree -> Branch("TheAoQ", &TheAoQ);
    //
    FrsTree -> Branch("Tof_wTref_S2_Cave", &Tof_wTref_S2_Cave);
    FrsTree -> Branch("Beta_S2_Cave", &Beta_S2_Cave);
    FrsTree -> Branch("Brho_S2_Cave", &Brho_S2_Cave);
    FrsTree -> Branch("AoQ_S2_Cave", &AoQ_S2_Cave);
    //
    FrsTree -> Branch("Tof_wTref_S2_S8", &Tof_wTref_S2_S8);
    FrsTree -> Branch("Beta_S2_S8", &Beta_S2_S8);
    FrsTree -> Branch("Brho_S2_S8", &Brho_S2_S8);
    FrsTree -> Branch("AoQ_S2_S8", &AoQ_S2_S8);
    //
    FrsTree -> Branch("Tof_wTref_S8_Cave", &Tof_wTref_S8_Cave);
    FrsTree -> Branch("Beta_S8_Cave", &Beta_S8_Cave);
    FrsTree -> Branch("Brho_S8_Cave", &Brho_S8_Cave);
    FrsTree -> Branch("AoQ_S8_Cave", &AoQ_S8_Cave);
    //
    return kSUCCESS;
}

void R3BSofFrsFillTree::Reset_Histo()
{
    LOG(INFO) << "R3BSofFrsFillTree::Reset_Histo";
}

void R3BSofFrsFillTree::Exec(Option_t* option)
{
    fNEvents += 1;
    //
    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(FATAL) << "R3BSofFrsFillTree::Exec FairRootManager not found";

    tpat = header->GetTpat();
    trigger = header->GetTrigger();
    if((header->GetTpat() & 1) == 0 && header->GetTpat() != 0) return; // Phisics events should not be registered as tpat==0 but there are some..

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

    if(MusicE< 0 && TwimE<0) return; // End this event to reduce output root file size
    
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
            if (!hitmapped) continue;
            iDet = hitmapped->GetDetector() - 1;
            iCh = hitmapped->GetPmt() - 1;
            multMapSci[iDet * fNbChannels + iCh]++;
            //fh1_finetime[iDet * fNbChannels + iCh]->Fill(hitmapped->GetTimeFine());
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
                if (!hittcal) continue;
                iDet = hittcal->GetDetector() - 1;
                iCh = hittcal->GetPmt() - 1;
                iRawTimeNs[iDet * fNbChannels + iCh] = hittcal->GetRawTimeNs();
            }
        }

        // --- ------------------------------ --- //
        // --- loop over sci single tcal data --- //
        // --- ------------------------------ --- //
	Int_t d,t;
        if (fSingleTcalItemsSci)
        {
            nHits = fSingleTcalItemsSci->GetEntriesFast();
            for (Int_t ihit = 0; ihit < nHits; ihit++)
            {
                R3BSofSciSingleTcalData* hitsingletcal = (R3BSofSciSingleTcalData*)fSingleTcalItemsSci->At(ihit);
                if (!hitsingletcal)   continue;
		d = hitsingletcal->GetDetector()-1;
		if (d == fIdS2-1)         xs2 = hitsingletcal->GetRawPosNs() * fS2SciCoef1 + fS2SciCoef0;
		if (d > 0) xpos[d-1] = hitsingletcal->GetRawPosNs();
		if(fIdS2>0&&fIdS8>0){
		  if(d==2) Tof_wTref_S2_S8 = hitsingletcal->GetRawTofNs_FromS2();
		  if(d==3){
		    Tof_wTref_S2_Cave = hitsingletcal->GetRawTofNs_FromS2(); // same as toff
		    Tof_wTref_S8_Cave = hitsingletcal->GetRawTofNs_FromS8();
		  }
		}
            }// end of loop over the SingleTcalItems

	    if (MusicE > 0){ // && xs2!=-5000.){ //&& toff!=-5000.) {
	      if(Tof_wTref_S2_Cave > 0.){
		Beta_S2_Cave = 462.837731 / (Tof_wTref_S2_Cave -1318.258541); // ToFCalib
		Gamma_S2_Cave = 1. / (TMath::Sqrt(1. - (Beta_S2_Cave) * (Beta_S2_Cave)));
		Brho_S2_Cave = fBrho0 * (1 + xs2 / 726.); //+mwpc0x/10./2000);
		AoQ_S2_Cave = Brho_S2_Cave / (3.10716 * Gamma_S2_Cave * Beta_S2_Cave);
	      }
	      //
	      if(Tof_wTref_S2_S8 > 0.){
		Beta_S2_S8 = 279.088230 / (Tof_wTref_S2_S8 -694.519095); // ToFCalib
		Gamma_S2_S8 = 1. / (TMath::Sqrt(1. - (Beta_S2_S8) * (Beta_S2_S8)));
		Brho_S2_S8 = fBrho0 * (1 + xs2 / 726.); //+mwpc0x/10./2000);
		AoQ_S2_S8 = Brho_S2_S8 / (3.10716 * Gamma_S2_S8 * Beta_S2_S8);
	      }
	      //
	      if(Tof_wTref_S8_Cave > 0.){
		Beta_S8_Cave = 183.845298 / (Tof_wTref_S8_Cave -623.62812); // ToFCalib
		Gamma_S8_Cave = 1. / (TMath::Sqrt(1. - (Beta_S8_Cave) * (Beta_S8_Cave)));
		Brho_S8_Cave = fBrho0 * (1 + (xs2 > -5000. ? xs2 : 0.) / 726.); //+mwpc0x/10./2000);
		AoQ_S8_Cave = Brho_S8_Cave / (3.10716 * Gamma_S8_Cave * Beta_S8_Cave);
	      }
	      //
	      if(Beta_S2_Cave > 0.){
		TheBeta = Beta_S2_Cave; TheGamma = Gamma_S2_Cave; TheBrho = Brho_S2_Cave; TheAoQ = AoQ_S2_Cave;
	      }else if(Beta_S2_S8 > 0.){
		TheBeta = Beta_S2_S8; TheGamma = Gamma_S2_S8; TheBrho = Brho_S2_S8; TheAoQ = AoQ_S2_S8;
	      }else if(Beta_S8_Cave > 0.){
		TheBeta = Beta_S8_Cave; TheGamma = Gamma_S8_Cave; TheBrho = Brho_S8_Cave; TheAoQ = AoQ_S8_Cave;
	      }
	      if(TheBeta > 0.){
		MusicZ = fZ0 + fZ1 * TMath::Sqrt(MusicE) * TheBeta
		  + fZ2 * MusicE * TheBeta * TheBeta;
		TwimZ = fTwimZ0 + fTwimZ1 * TMath::Sqrt(TwimE) * TheBeta
		  + fTwimZ2 * TwimE * TheBeta * TheBeta;
	      }
	      //
	    } //if MusicE>0
	}
    }

    //
    FrsTree->Fill();
}

// -----   Public method Reset   ------------------------------------------------
void R3BSofFrsFillTree::Reset()
{
}

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
  //Init branch values
  MusicZ = -5000., MusicE = -5000.;
  MusicDT = -5000.;
  TwimE = -5000., TwimZ = -5000.;
  xs2 = -5000.;
  xpos[0] = -5000.;  xpos[1] = -5000.;  xpos[2] = -5000.;
  //toff = -5000.;
  Tof_wTref_S2_Cave = -5000., Beta_S2_Cave = -5000., Gamma_S2_Cave = -5000., Brho_S2_Cave = -5000.;
  Tof_wTref_S2_S8 = -5000., Beta_S2_S8 = -5000., Gamma_S2_S8 = -5000., Brho_S2_S8 = -5000.;
  Tof_wTref_S8_Cave = -5000., Beta_S8_Cave = -5000., Gamma_S8_Cave = -5000., Brho_S8_Cave = -5000.;
  AoQ_S2_Cave = -5000., AoQ_S2_S8 = -5000., AoQ_S8_Cave = -5000.;
  TheBeta = -5000., TheGamma = -5000., TheBrho = -5000., TheAoQ = -5000.;
  tpat = 0, trigger =0;
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

void R3BSofFrsFillTree::FinishTask()
{
  FrsTree->Write();
}

ClassImp(R3BSofFrsFillTree)
