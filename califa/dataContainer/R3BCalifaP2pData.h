/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019 Members of R3B Collaboration                          *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

#ifndef R3BCALIFAP2PDATA_H
#define R3BCALIFAP2PDATA_H

#include "TObject.h"

#include "FairMultiLinkedData.h"
#include "R3BCalifaCrystalCalData.h"

class R3BCalifaP2pData : public FairMultiLinkedData
{

  public:
    /** Default constructor **/
    R3BCalifaP2pData();

    /** Constructor with arguments
     *@param fEnergy1 Energy of proton/neutron hit candidate 1
     *@param fEnergy2 Energy of proton/neutron hit candidate 1
     *@param fEnergySum The sum energy of the two hit candidates
     *@param fPhi1 The phi angle of proton/neutron hit candidate 1
     *@param fPhi2 The phi angle of proton/neutron hit candidate 2
     *@param fPhiSum The positive phi angle between the two candidates
     *@param fTheta1 The theta angle of proton/neutron hit candidate 1
     *@param fTheta2 The theta angle of proton/neutron hit candidate 2
     *@param fThetaSum The positive theta angle between the two candidates. Only correct, when phi = 180 degree!
     
     
     **/
    R3BCalifaP2pData(Double_t ene1, 
                    Double_t ene2, 
                    Double_t eneSum, 
                    Double_t phi1, 
                    Double_t phi2, 
                    Double_t phiSum, 
                    Double_t theta1, 
                    Double_t theta2, 
                    Double_t thetaSum, 
                    ULong64_t time);
          
    /** Destructor **/
    virtual ~R3BCalifaP2pData() {}

    /** Accessors **/
    Double_t GetEnergy1() const { return fEnergy1; }
    Double_t GetEnergy2() const { return fEnergy2; }
    Double_t GetEnergySum() const { return fEnergySum; }
    Double_t GetTheta1() const { return fTheta1; }
    Double_t GetTheta2() const { return fTheta2; }
    Double_t GetThetaSum() const { return fThetaSum; }
    Double_t GetPhi1() const { return fPhi1; }
    Double_t GetPhi2() const { return fPhi2; }
    Double_t GetPhiSum() const { return fPhiSum; }
    ULong64_t GetTime() const { return fTime; }

    /** Modifiers **/
    void SetEnergy1(Double_t ene) { fEnergy1 = ene; }
    void SetEnergy2(Double_t ene) { fEnergy2 = ene; }
    void SetEnergySum(Double_t ene) { fEnergySum = ene; }
    void SetTheta1(Double_t theta) { fTheta1 = theta; }
    void SetTheta2(Double_t theta) { fTheta2 = theta; }
    void SetThetaSum(Double_t theta) { fThetaSum = theta; }
    void SetPhi1(Double_t phi) { fPhi1 = phi; }
    void SetPhi2(Double_t phi) { fPhi2 = phi; }
    void SetPhiSum(Double_t phi) { fPhiSum = phi; }
    void SetTime(ULong64_t time) { fTime = time; }

  protected:
    // Basic Hit information
    Double_t fEnergy1;        // Energy deposition particle 1
    Double_t fEnergy2;        // Energy deposition particle 2
    Double_t fEnergySum;      // Sum energy of the two particles

    Double_t fTheta1;         // Theta of particle 1
    Double_t fTheta2;         // Theta of particle 2
    Double_t fThetaSum;       // Theta between the two particles, if phi=180 degree!
    
    Double_t fPhi1;           // Phi of particle 1
    Double_t fPhi2;           // Phi of particle 2
    Double_t fPhiSum;         // Phi between the two particles
    
    ULong64_t fTime;          // WR time stamp
    

    ClassDef(R3BCalifaP2pData, 5)
};

#endif
