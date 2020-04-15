// ------------------------------------------------------------
// -----                  R3BSofSciOnlineSpectra           -----
// -----           Fill SOFIA online histograms           -----
// ------------------------------------------------------------

/*
 * This task should fill histograms with SOFIA online data
 */

#include "R3BSofSciOnlineSpectra.h"
#include "R3BEventHeader.h"
#include "R3BMusicCalData.h"
#include "R3BMusicHitData.h"
#include "R3BSofMwpcCalData.h"
#include "R3BSofMwpcHitData.h"
#include "R3BSofSciMappedData.h"
#include "R3BSofSciCalData.h"
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

R3BSofSciOnlineSpectra::R3BSofSciOnlineSpectra()
    : FairTask("SofSciOnlineSpectra", 1)
    , fMappedItemsSci(NULL)
    , fTcalItemsSci(NULL)
    , fMusHitItems(NULL)
    , fMusCalItems(NULL)
    , fCalItemsMwpc0(NULL)
    , fTofwHitData(NULL)
    , fNEvents(0)
{
}

R3BSofSciOnlineSpectra::R3BSofSciOnlineSpectra(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fMappedItemsSci(NULL)
    , fTcalItemsSci(NULL)
    , fMusHitItems(NULL)
    , fMusCalItems(NULL)
    , fCalItemsMwpc0(NULL)
    , fTofwHitData(NULL)
    , fNEvents(0)
{
}

R3BSofSciOnlineSpectra::~R3BSofSciOnlineSpectra()
{
    LOG(INFO) << "R3BSofSciOnlineSpectra::Delete instance";
    if (fMappedItemsSci)
        delete fMappedItemsSci;
    if (fTcalItemsSci)
        delete fTcalItemsSci;
    if (fMusHitItems)
        delete fMusHitItems;
    if (fMusCalItems)
        delete fMusCalItems;
    if (fCalItemsMwpc0)
        delete fCalItemsMwpc0;
    if (fTofwHitData)
        delete fTofwHitData;
}

InitStatus R3BSofSciOnlineSpectra::Init()
{

    LOG(INFO) << "R3BSofSciOnlineSpectra::Init ";

    // try to get a handle on the EventHeader. EventHeader may not be
    // present though and hence may be null. Take care when using.

    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(FATAL) << "R3BSofSciOnlineSpectra::Init FairRootManager not found";
    // header = (R3BEventHeader*)mgr->GetObject("R3BEventHeader");

    FairRunOnline* run = FairRunOnline::Instance();
    run->GetHttpServer()->Register("", this);

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

    // get access to hit data of the MUSIC
    fMusHitItems = (TClonesArray*)mgr->GetObject("MusicHitData");
    if (!fMusHitItems)
        LOG(WARNING) << "R3BSofSciOnlineSpectra: MusicHitData not found";

    // get access to cal data of the MUSIC
    fMusCalItems = (TClonesArray*)mgr->GetObject("MusicCalData");
    if (!fMusCalItems)
        LOG(WARNING) << "R3BSofSciOnlineSpectra: MusicCalData not found";

    // get access to cal data of the MWPC0
    //fCalItemsMwpc0 = (TClonesArray*)mgr->GetObject("Mwpc0CalData");
    //if (!fCalItemsMwpc0)
    //    LOG(WARNING) << "R3BSofSciOnlineSpectra: Mwpc0CalData not found";
    fCalItemsMwpc0 = (TClonesArray*)mgr->GetObject("Mwpc0HitData");
    if (!fCalItemsMwpc0)
        LOG(WARNING) << "R3BSofSciOnlineSpectra: Mwpc0HitData not found";


    // --- ------------------------------- --- //
    // --- Create histograms for detectors --- //
    // --- ------------------------------- --- //
    char Name1[255];
    char Name2[255];

    for (Int_t i = 0; i < NbDetectors; i++)
    {

        // === FINE TIME AND MULT === //
        sprintf(Name1, "SofSci%i_MultAndFt", i + 1);
        cSciMult[i] = new TCanvas(Name1, Name1, 10, 10, 800, 700);
        cSciMult[i]->Divide(2, 2);
        for (Int_t j = 0; j < NbChannels; j++)
        {
            sprintf(Name1, "SofSci%i_FineTimeCh%i", i + 1, j + 1);
            fh1_finetime[i * NbChannels + j] = new TH1I(Name1, Name1, 1000, 0, 1000);
            cSciMult[i]->cd(j + 1);
            fh1_finetime[i * NbChannels + j]->Draw("");
        }
        sprintf(Name1, "SofSci%i_MultPerChannel", i + 1);
        fh2_mult[i] = new TH2I(Name1, Name1, NbChannels, 0.5, NbChannels + 0.5, 20, -0.5, 19.5);
        fh2_mult[i]->GetXaxis()->SetTitle("channel: 1=PMT R,    2=PMT L,    3=COMMON REF");
        fh2_mult[i]->GetYaxis()->SetTitle("multiplicity per channel");
        cSciMult[i]->cd(4);
        fh2_mult[i]->Draw("COL");

        // === RAW POSITION === //
        sprintf(Name1, "SofSci%i_RawPos", i + 1);
        cSciRawPos[i] = new TCanvas(Name1, Name1, 10, 10, 500, 500);
        sprintf(Name1, "SofSci%i_RawPosAtTcal_Mult1", i + 1);
        fh1_RawPos_AtTcalMult1[i] = new TH1F(Name1, Name1, 10000, -10, 10);
        fh1_RawPos_AtTcalMult1[i]->GetXaxis()->SetTitle("Raw position [ns with one bin/ps]");
        fh1_RawPos_AtTcalMult1[i]->GetYaxis()->SetTitle("Counts per bin");
        fh1_RawPos_AtTcalMult1[i]->GetXaxis()->CenterTitle(true);
        fh1_RawPos_AtTcalMult1[i]->GetYaxis()->CenterTitle(true);
        fh1_RawPos_AtTcalMult1[i]->GetXaxis()->SetLabelSize(0.045);
        fh1_RawPos_AtTcalMult1[i]->GetXaxis()->SetTitleSize(0.045);
        fh1_RawPos_AtTcalMult1[i]->GetYaxis()->SetLabelSize(0.045);
        fh1_RawPos_AtTcalMult1[i]->GetYaxis()->SetTitleSize(0.045);
        cSciRawPos[i]->cd();
        fh1_RawPos_AtTcalMult1[i]->Draw("");
    }
  Int_t nTof=0;
  for (Int_t dstart = 0; dstart < NbDetectors-1 ; dstart++)
  {
    for(Int_t dstop = dstart+1; dstop < NbDetectors ; dstop++)
    {
      //1D - Raw time-of-flight from tcal data at mult=1
      sprintf(Name1,"RawTof_Sci%02d_to_Sci%02d",dstart+1,dstop+1);
      sprintf(Name2,"Raw time-of-flight from Sci%02d to Sci%02d",dstart+1, dstop+1);
      cSciRawTof[nTof] = new TCanvas(Name1,Name2,10,10,800,800);
      cSciRawTof[nTof]->Divide(1,2);
      sprintf(Name1,"RawTofNs_Sci%02d_to_Sci%02d",dstart+1,dstop+1);
      fh1_RawTof_AtTcalMult1[nTof] = new TH1D(Name1,Name1,1000000,-100000,100000);
      fh1_RawTof_AtTcalMult1[nTof]->GetXaxis()->SetTitle("Raw Tof [ns]");
      fh1_RawTof_AtTcalMult1[nTof]->GetYaxis()->SetTitle("Counts per bin"); 
      sprintf(Name1,"RawTofNs_wTref_Sci%02d_to_Sci%02d",dstart+1,dstop+1);
      fh1_RawTof_AtTcalMult1_wTref[nTof] = new TH1D(Name1,Name1,100000,1000,2000);
      fh1_RawTof_AtTcalMult1_wTref[nTof]->GetXaxis()->SetTitle("Raw Tof [ns]");
      fh1_RawTof_AtTcalMult1_wTref[nTof]->GetYaxis()->SetTitle("Counts per bin"); 
      cSciRawTof[nTof]->cd(1);
      fh1_RawTof_AtTcalMult1[nTof]->Draw("");
      cSciRawTof[nTof]->cd(2);
      fh1_RawTof_AtTcalMult1_wTref[nTof]->Draw("");
     
      // Raw Tof versus position in start and stop
      sprintf(Name1,"RawTof_Sci%02d_Sci%02d_vs_RawPos",dstart+1,dstop+1);
      sprintf(Name2,"Raw time-of-flight from Sci%02d to Sci%02d versus RawPos",dstart+1, dstop+1);
      cSciRawTofvsRawPos[nTof] = new TCanvas(Name1,Name2,10,10,800,800);
      cSciRawTofvsRawPos[nTof]->Divide(1,2);
      sprintf(Name1,"RawTofNs_Sci%02d_Sci%02d_vs_RawPos%02d",dstart+1,dstop+1,dstart+1);
      fh2_RawTof_vs_RawPosStart_AtTcalMult1[nTof] = new TH2D(Name1,Name1,1000,-10,10,2000,-10000,1000);
      fh2_RawTof_vs_RawPosStart_AtTcalMult1[nTof]->GetXaxis()->SetTitle("Raw Pos [ns]");
      fh2_RawTof_vs_RawPosStart_AtTcalMult1[nTof]->GetYaxis()->SetTitle("Raw Tof [ns]"); 
      sprintf(Name1,"RawTofNs_Sci%02d_Sci%02d_vs_RawPos%02d",dstart+1,dstop+1,dstop+1);
      fh2_RawTof_vs_RawPosStop_AtTcalMult1[nTof] = new TH2D(Name1,Name1,1000,-10,10,2000,-10000,1000);
      fh2_RawTof_vs_RawPosStop_AtTcalMult1[nTof]->GetXaxis()->SetTitle("Raw Pos [ns]");
      fh2_RawTof_vs_RawPosStop_AtTcalMult1[nTof]->GetYaxis()->SetTitle("Raw Tof [ns]"); 
      cSciRawTofvsRawPos[nTof]->cd(1);
      fh2_RawTof_vs_RawPosStart_AtTcalMult1[nTof]->Draw("col");
      cSciRawTofvsRawPos[nTof]->cd(2);
      fh2_RawTof_vs_RawPosStop_AtTcalMult1[nTof]->Draw("col");
      nTof++;
    }  
  }

    // Music Hit data vs SCI-RawPos
    TCanvas* cMusicZvsRawPos =
        new TCanvas("Musicchargez_vs_RawPosAtTcal_Mult1", "Music charge Z vs RawPosAtTcal_Mult1", 10, 10, 800, 700);
    fh2_MusZvsRawPos =
        new TH2F("fh2_Musicchargez_vs_RawPos", "Music charge Z vs RawPosAtTcal_Mult1", 10000, -7, 7, 1200, 6, 40);
    fh2_MusZvsRawPos->GetXaxis()->SetTitle("Raw position [ns with one bin/ps]");
    fh2_MusZvsRawPos->GetYaxis()->SetTitle("Charge (Z)");
    fh2_MusZvsRawPos->GetYaxis()->SetTitleOffset(1.1);
    fh2_MusZvsRawPos->GetXaxis()->CenterTitle(true);
    fh2_MusZvsRawPos->GetYaxis()->CenterTitle(true);
    fh2_MusZvsRawPos->GetXaxis()->SetLabelSize(0.045);
    fh2_MusZvsRawPos->GetXaxis()->SetTitleSize(0.045);
    fh2_MusZvsRawPos->GetYaxis()->SetLabelSize(0.045);
    fh2_MusZvsRawPos->GetYaxis()->SetTitleSize(0.045);
    fh2_MusZvsRawPos->Draw("colz");


    // Music Hit data vs SCI-RawTof
    TCanvas* cMusicZvsRawTof =
        new TCanvas("Music_chargeZ_vs_RawTofAtTcal_Mult1", "Music charge Z vs RawTofAtTcal_Mult1", 10, 10, 800, 700);
    fh2_MusZvsRawTof =
      //new TH2F("fh2_Music_chargez_vs_RawTof", "Music charge Z vs RawTofAtTcal_Mult1 (Cave-C)", 4000, -3600, -2900, 1200, 6, 40);
        new TH2F("fh2_Music_chargez_vs_RawTof", "Music charge Z vs RawTofAtTcal_Mult1 (Cave-C)", 4000, -10000, 10000, 1200, 6, 40);
    fh2_MusZvsRawTof->GetXaxis()->SetTitle("Raw-ToF-Cave-C [ns]");
    fh2_MusZvsRawTof->GetYaxis()->SetTitle("Charge (Z)");
    fh2_MusZvsRawTof->GetYaxis()->SetTitleOffset(1.1);
    fh2_MusZvsRawTof->GetXaxis()->CenterTitle(true);
    fh2_MusZvsRawTof->GetYaxis()->CenterTitle(true);
    fh2_MusZvsRawTof->GetXaxis()->SetLabelSize(0.045);
    fh2_MusZvsRawTof->GetXaxis()->SetTitleSize(0.045);
    fh2_MusZvsRawTof->GetYaxis()->SetLabelSize(0.045);
    fh2_MusZvsRawTof->GetYaxis()->SetTitleSize(0.045);
    fh2_MusZvsRawTof->Draw("colz");

    // Hit data, Aq_vs_q
    TCanvas* cAqvsq = new TCanvas("FRSv_Aq_vs_q", "A/q_vs_q 2D info", 10, 10, 800, 700);
    fh2_Aqvsq = new TH2F("fh2v_Aq_vs_q_frs", "FRS: A/q vs q", 3000, 1., 3, 1300, 8, 39.5);
    fh2_Aqvsq->GetXaxis()->SetTitle("A/q");
    fh2_Aqvsq->GetYaxis()->SetTitle("Z [Charge units]");
    fh2_Aqvsq->GetYaxis()->SetTitleOffset(1.1);
    fh2_Aqvsq->GetXaxis()->CenterTitle(true);
    fh2_Aqvsq->GetYaxis()->CenterTitle(true);
    fh2_Aqvsq->GetXaxis()->SetLabelSize(0.045);
    fh2_Aqvsq->GetXaxis()->SetTitleSize(0.045);
    fh2_Aqvsq->GetYaxis()->SetLabelSize(0.045);
    fh2_Aqvsq->GetYaxis()->SetTitleSize(0.045);
    fh2_Aqvsq->Draw("colz");


    // Music Hit data vs SCI-RawTof-S8
    TCanvas* cMusicZvsRawTofS8 =
        new TCanvas("Music_chargeZ_vs_RawTofAtTcal_Mult1S8", "Music charge Z vs RawTofAtTcal_Mult1 S8", 10, 10, 800, 700);
    fh2_MusZvsRawTofS8 =
        new TH2F("fh2_Music_chargez_vs_RawTof_S8", "Music charge Z vs RawTofAtTcal_Mult1 (S8)", 20000, -1000, 1000, 1200, 6, 40);
    fh2_MusZvsRawTofS8->GetXaxis()->SetTitle("Raw-ToF-S8 [ns]");
    fh2_MusZvsRawTofS8->GetYaxis()->SetTitle("Charge (Z)");
    fh2_MusZvsRawTofS8->GetYaxis()->SetTitleOffset(1.1);
    fh2_MusZvsRawTofS8->GetXaxis()->CenterTitle(true);
    fh2_MusZvsRawTofS8->GetYaxis()->CenterTitle(true);
    fh2_MusZvsRawTofS8->GetXaxis()->SetLabelSize(0.045);
    fh2_MusZvsRawTofS8->GetXaxis()->SetTitleSize(0.045);
    fh2_MusZvsRawTofS8->GetYaxis()->SetLabelSize(0.045);
    fh2_MusZvsRawTofS8->GetYaxis()->SetTitleSize(0.045);
    fh2_MusZvsRawTofS8->Draw("colz");

    // Music Hit data vs SCI-RawPos
    TCanvas* cMusicDTvsRawPos =
        new TCanvas("MusicDT_vs_RawPosAtTcal_Mult1", "Music FT vs RawPosAtTcal_Mult1", 10, 10, 800, 700);
    fh2_MusDTvsRawPos = new TH2F("fh2_MusicDT_vs_RawPos", "Music DT vs RawPosAtTcal_Mult1", 10000, -7, 7, 800, -20, 20);
    fh2_MusDTvsRawPos->GetXaxis()->SetTitle("Sci Raw position [ns with one bin/ps]");
    fh2_MusDTvsRawPos->GetYaxis()->SetTitle("Drift Time (mm)");
    fh2_MusDTvsRawPos->GetYaxis()->SetTitleOffset(1.1);
    fh2_MusDTvsRawPos->GetXaxis()->CenterTitle(true);
    fh2_MusDTvsRawPos->GetYaxis()->CenterTitle(true);
    fh2_MusDTvsRawPos->GetXaxis()->SetLabelSize(0.045);
    fh2_MusDTvsRawPos->GetXaxis()->SetTitleSize(0.045);
    fh2_MusDTvsRawPos->GetYaxis()->SetLabelSize(0.045);
    fh2_MusDTvsRawPos->GetYaxis()->SetTitleSize(0.045);
    fh2_MusDTvsRawPos->Draw("col");

    // Mwpc0 cal data vs SCI-RawPos
    TCanvas* cMwpc0vsRawPos =
        new TCanvas("Mwpc0_vs_RawPosAtTcal_Mult1", "Mwpc0-X vs RawPosAtTcal_Mult1", 10, 10, 800, 700);
    fh2_Mwpc0vsRawPos = new TH2F("fh2_Mwpc_vs_RawPos", "Mwpc0-X vs RawPosAtTcal_Mult1", 10000, -7, 7, 400, -100, 100);
    fh2_Mwpc0vsRawPos->GetXaxis()->SetTitle("Raw position [ns with one bin/ps]");
    fh2_Mwpc0vsRawPos->GetYaxis()->SetTitle("Mwpc0-X [pads]");
    fh2_Mwpc0vsRawPos->GetYaxis()->SetTitleOffset(1.1);
    fh2_Mwpc0vsRawPos->GetXaxis()->CenterTitle(true);
    fh2_Mwpc0vsRawPos->GetYaxis()->CenterTitle(true);
    fh2_Mwpc0vsRawPos->GetXaxis()->SetLabelSize(0.045);
    fh2_Mwpc0vsRawPos->GetXaxis()->SetTitleSize(0.045);
    fh2_Mwpc0vsRawPos->GetYaxis()->SetLabelSize(0.045);
    fh2_Mwpc0vsRawPos->GetYaxis()->SetTitleSize(0.045);
    fh2_Mwpc0vsRawPos->Draw("col");

    // --- --------------- --- //
    // --- MAIN FOLDER-Sci --- //
    // --- --------------- --- //
    TFolder* mainfolSci = new TFolder("SOFSCI", "SOFSCI info");
    for (Int_t i = 0; i < NbDetectors; i++)
    {
        mainfolSci->Add(cSciMult[i]);
        mainfolSci->Add(cSciRawPos[i]);
    }
    for(Int_t i = 0; i < NbTof ; i++)
    {
	mainfolSci->Add(cSciRawTof[i]);
    }
    mainfolSci->Add(cMusicDTvsRawPos);
    mainfolSci->Add(cMwpc0vsRawPos);
    run->AddObject(mainfolSci);

    TFolder* mainfolID = new TFolder("IncomingID", "IncomingID info");
    mainfolID->Add(cMusicZvsRawPos);
    mainfolID->Add(cMusicZvsRawTof);
    mainfolID->Add(cMusicZvsRawTofS8);
    mainfolID->Add(cAqvsq);
    run->AddObject(mainfolID);

    // Register command to reset histograms
    run->GetHttpServer()->RegisterCommand("Reset_SOFSCI_HIST", Form("/Objects/%s/->Reset_Histo()", GetName()));



    // OUTPUT DATA
    fTofwHitData = new TClonesArray("R3BSofSciCalData", 5);
 
        //mgr->Register("SofSciHitData", "Tof FRS Hit", fTofwHitData, kTRUE);
    mgr->Register("SofSciHitData", "Tof FRS Hit", fTofwHitData, kFALSE);
    
    return kSUCCESS;
}

void R3BSofSciOnlineSpectra::Reset_Histo()
{
    LOG(INFO) << "R3BSofSciOnlineSpectra::Reset_Histo";
    for (Int_t i = 0; i < NbDetectors; i++)
    {
        // === MULT AND FINE TIME === //
        fh2_mult[i]->Reset();
        for (Int_t j = 0; j < NbChannels; j++)
        {
            fh1_finetime[i * NbChannels + j]->Reset();
        }
        // === RAW POSITION === //
        fh1_RawPos_AtTcalMult1[i]->Reset();
    }
    for(Int_t i = 0; i<NbTof; i++)
    {
      // === RAW TIME_OF_FLIGHT === //
      fh1_RawTof_AtTcalMult1[i]->Reset();
      fh1_RawTof_AtTcalMult1_wTref[i]->Reset();
      fh2_RawTof_vs_RawPosStart_AtTcalMult1[i]->Reset();
      fh2_RawTof_vs_RawPosStop_AtTcalMult1[i]->Reset();
    }
    fh2_MusZvsRawPos->Reset();
    fh2_MusZvsRawTof->Reset();
    fh2_MusZvsRawTofS8->Reset();
    fh2_Aqvsq->Reset();
    fh2_MusDTvsRawPos->Reset();
    fh2_Mwpc0vsRawPos->Reset();
}

void R3BSofSciOnlineSpectra::Exec(Option_t* option)
{
    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(FATAL) << "R3BSofSciOnlineSpectra::Exec FairRootManager not found";

    Int_t nHits;
    UShort_t iDet; // 0-bsed
    UShort_t iCh;  // 0-based
    Double_t iRawTimeNs[NbDetectors * NbChannels];
    UShort_t mult[NbDetectors * NbChannels];

    // --- -------------- --- //
    // --- initialisation --- //
    // --- -------------- --- //
    for (UShort_t i = 0; i < NbDetectors; i++)
    {
        for (UShort_t j = 0; j < NbChannels; j++)
        {
            mult[i * NbChannels + j] = 0;
        }
    }

    // MUSIC Hit data
    Double_t MusicZ = 0.;
    Double_t MusicDT = -1000000.;
    if (fMusHitItems && fMusHitItems->GetEntriesFast() > 0)
    {
        nHits = fMusHitItems->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BMusicHitData* hit = (R3BMusicHitData*)fMusHitItems->At(ihit);
            if (!hit)
                continue;
            MusicZ = hit->GetZcharge();
        }
    }
    // MUSIC Cal data
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
    }

    if (fMappedItemsSci && fMappedItemsSci->GetEntriesFast() && fTcalItemsSci && fTcalItemsSci->GetEntriesFast())
    {

        // --- --------------------- --- //
        // --- loop over mapped data --- //
        // --- --------------------- --- //
        nHits = fMappedItemsSci->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BSofSciMappedData* hitmapped = (R3BSofSciMappedData*)fMappedItemsSci->At(ihit);
            if (!hitmapped)
                continue;
            iDet = hitmapped->GetDetector() - 1;
            iCh = hitmapped->GetPmt() - 1;
            mult[iDet * NbChannels + iCh]++;
            fh1_finetime[iDet * NbChannels + iCh]->Fill(hitmapped->GetTimeFine());
        }

        // --- ------------------- --- //
        // --- loop over tcal data --- //
        // --- ----------fTofwHitData--------- --- //
        nHits = fTcalItemsSci->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BSofSciTcalData* hittcal = (R3BSofSciTcalData*)fTcalItemsSci->At(ihit);
            if (!hittcal)
                continue;
            iDet = hittcal->GetDetector() - 1;
            iCh = hittcal->GetPmt() - 1;
            iRawTimeNs[iDet * NbChannels + iCh] = hittcal->GetRawTimeNs();
            ;
        }

        // Get cal data MWPC0
        Double_t mwpc0x = -1000., qmax = -100.;
        if (fCalItemsMwpc0 && fCalItemsMwpc0->GetEntriesFast() > 0)
        {
            nHits = fCalItemsMwpc0->GetEntriesFast();
            for (Int_t ihit = 0; ihit < nHits; ihit++)
            {
                //R3BSofMwpcCalData* hit = (R3BSofMwpcCalData*)fCalItemsMwpc0->At(ihit);
                R3BSofMwpcHitData* hit = (R3BSofMwpcHitData*)fCalItemsMwpc0->At(ihit);
                if (!hit)
                    continue;
		mwpc0x = hit->GetX();
                //if (hit->GetQ() > qmax && hit->GetPlane() == 1)
                //{
                //    mwpc0x = hit->GetX();
                //    qmax = hit->GetQ();
                //}
            }
        }

        // --- ----------------------------------------- --- //
        // --- filling some histogramms outside the loop --- //
        // --- ----------------------------------------- --- //
        Double_t possci = 0.;
        Double_t xs2=-10000.;  
        Double_t xcaveC= -10000;	
        for (UShort_t i = 0; i < NbDetectors; i++)
        {
            for (UShort_t j = 0; j < NbChannels; j++)
            {
                fh2_mult[i]->Fill(j + 1, mult[i * NbChannels + j]);
            }
            if ((mult[i * NbChannels] > 0) && (mult[i * NbChannels + 1] > 0))
            {
                // x position increases from left to right : TrawRIGHT - TrawLEFT
                possci = iRawTimeNs[i * NbChannels] - iRawTimeNs[i * NbChannels + 1];
		if(i==3) xcaveC = possci;
                fh1_RawPos_AtTcalMult1[i]->Fill(possci);
		Double_t slope_calib = -5.8;
                if(i==1)xs2=possci*slope_calib;
                if (MusicZ > 0. && i==3)
                    fh2_MusZvsRawPos->Fill(possci, MusicZ);

                if (MusicDT != -1000000. && i==3)
                    fh2_MusDTvsRawPos->Fill(possci, MusicZ);

                if (mwpc0x != -1000. && possci > -10. && possci < 10. && i==3)
                {
                    fh2_Mwpc0vsRawPos->Fill(possci, mwpc0x);
                }
            }
        }
	Int_t indexTof=0;
	Double_t iRawTof;
	for(UShort_t dstart=0; dstart<NbDetectors-1; dstart++)
	{
	  for(UShort_t dstop=dstart+1; dstop<NbDetectors; dstop++)
	  {
	    if( (mult[dstart*NbChannels] == 1) && (mult[dstart*NbChannels+1] == 1) &&
		(mult[dstop*NbChannels] == 1) && (mult[dstop*NbChannels+1] == 1) )
	    {
	      iRawTof = 0.5*(iRawTimeNs[dstop*NbChannels]+iRawTimeNs[dstop*NbChannels+1]) -
			0.5*(iRawTimeNs[dstart*NbChannels]+iRawTimeNs[dstart*NbChannels+1]) + 
			iRawTimeNs[dstart*NbChannels+2] - iRawTimeNs[dstop*NbChannels+2];
	      fh1_RawTof_AtTcalMult1_wTref[indexTof]->Fill(iRawTof); 
	      fh2_RawTof_vs_RawPosStart_AtTcalMult1[indexTof]->Fill(iRawTimeNs[dstart*NbChannels]-iRawTimeNs[dstart*NbChannels+1],iRawTof);
	      fh2_RawTof_vs_RawPosStop_AtTcalMult1[indexTof]->Fill(iRawTimeNs[dstop*NbChannels]-iRawTimeNs[dstop*NbChannels+1],iRawTof);
	      fh1_RawTof_AtTcalMult1[indexTof]->Fill( 
		  0.5*(iRawTimeNs[dstop*NbChannels]+iRawTimeNs[dstop*NbChannels+1]) - 
		  0.5*(iRawTimeNs[dstart*NbChannels]+iRawTimeNs[dstart*NbChannels+1])); 

              if(indexTof==4){ 
              double toff=0.5*(iRawTimeNs[dstop*NbChannels]+iRawTimeNs[dstop*NbChannels+1]) - 
		   0.5*(iRawTimeNs[dstart*NbChannels]+iRawTimeNs[dstart*NbChannels+1]);
		      
	      AddHitData(2, iRawTimeNs[dstop*NbChannels]-iRawTimeNs[dstop*NbChannels+1], toff);
	        //AddHitData(2, iRawTimeNs[dstop*NbChannels]-iRawTimeNs[dstop*NbChannels+1], iRawTof);
     
              fh2_MusZvsRawTofS8->Fill(0.5*(iRawTimeNs[dstop*NbChannels]+iRawTimeNs[dstop*NbChannels+1]) - 
		   0.5*(iRawTimeNs[dstart*NbChannels]+iRawTimeNs[dstart*NbChannels+1]), MusicZ);

                if(xs2>-900.)AddHitData(1, xs2, -50.);
              }
              if(indexTof==4){
              double toff=0.5*(iRawTimeNs[dstop*NbChannels]+iRawTimeNs[dstop*NbChannels+1]) - 
		   0.5*(iRawTimeNs[dstart*NbChannels]+iRawTimeNs[dstart*NbChannels+1]);
              fh2_MusZvsRawTof->Fill(toff, MusicZ);

	      //              double Beta_S2_Cave = 15424.3/(toff+675.+3224.)/29.9999; // this is for runs up to 303 and earlier
              //double Beta_S2_Cave = 15424.3/(toff+665.-3946.)/29.9999;//239
	      //double Beta_S2_Cave = 15424.3/(toff+675.+13443.)/29.9999;//until morning 24 Feb
	      //double Beta_S2_Cave = 15424.3/(toff+675.+12433.)/29.9999;// After run313
	     // double Beta_S2_Cave = 15424.3/(toff+675.+9778.)/29.9999;// After run 316
	      //	  double Beta_S2_Cave = 15424.3/(toff+675.+9779)/29.9999;// After run 336

	      //	  double Beta_S2_Cave = 15424.3/(toff+698.83+129.4)/29.9999;// for 38Ca runs 263-265
	      double Beta_S2_Cave = 15424.3/(toff+698.83-9011.52)/29.9999;// for 38Ca runs 263-265
	      double Gamma_S2_Cave = 1. / (TMath::Sqrt(1. - (Beta_S2_Cave) * (Beta_S2_Cave)));
	      //	      double Brho_S2_Cave = 9.048*(1+xs2/726);//+mwpc0x/10./2000);
	      double Brho_S2_Cave = 6.7422063*(1+xs2/667./*726.*/); // effective Brho for runs 263-265
              fh2_Aqvsq->Fill(Brho_S2_Cave/ (3.10716 * Gamma_S2_Cave * Beta_S2_Cave), MusicZ);
              }
	    }
	    indexTof++;
	  }
	}
    }
    fNEvents += 1;
}


// -----   Public method Reset   ------------------------------------------------
void R3BSofSciOnlineSpectra::Reset()
{
    LOG(DEBUG) << "Clearing TofWHitData Structure";
    if (fTofwHitData)
        fTofwHitData->Clear();
}

// -----   Private method AddHitData  --------------------------------------------
R3BSofSciHitData* R3BSofSciOnlineSpectra::AddHitData(Int_t paddle, Double_t X, Double_t time)
{
    // It fills the R3BSofTofwHitData
    TClonesArray& clref = *fTofwHitData;
    Int_t size = clref.GetEntriesFast();
    return new (clref[size]) R3BSofSciHitData(paddle, X, time);
}

void R3BSofSciOnlineSpectra::FinishEvent()
{
    /*if (fMappedItemsSci)
    {
        fMappedItemsSci->Clear();
    }*/
    if (fTcalItemsSci)
    {
        fTcalItemsSci->Clear();
    }
    /*if (fMusHitItems)
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
    }*/
}

void R3BSofSciOnlineSpectra::FinishTask()
{

    if (fMappedItemsSci)
    {
        for (UShort_t i = 0; i < NbDetectors; i++)
        {
            fh2_mult[i]->Write();
            cSciMult[i]->Write();
            for (UShort_t j = 0; j < NbChannels; j++)
            {
                fh1_finetime[i * NbChannels + j]->Write();
            }
        }
    }

    if (fTcalItemsSci)
    {
        for (UShort_t i = 0; i < NbDetectors; i++)
        {
            fh1_RawPos_AtTcalMult1[i]->Write();
            cSciRawPos[i]->Write();
        }
	for (UShort_t i = 0; i<NbTof ; i++)
	{
	    fh1_RawTof_AtTcalMult1[i]->Write();
	    fh1_RawTof_AtTcalMult1_wTref[i]->Write();
	    cSciRawTof[i]->Write();
	    fh2_RawTof_vs_RawPosStart_AtTcalMult1[i]->Write();
	    fh2_RawTof_vs_RawPosStop_AtTcalMult1[i]->Write();
	    cSciRawTofvsRawPos[i];
	}
        if (fMusHitItems)
        {
            fh2_MusZvsRawPos->Write();
        }
        if (fMusCalItems){
           fh2_MusZvsRawTof->Write();
           fh2_MusDTvsRawPos->Write();
           fh2_Aqvsq->Write();
        }
        if (fCalItemsMwpc0)
            fh2_Mwpc0vsRawPos->Write();
    }
}

ClassImp(R3BSofSciOnlineSpectra)
