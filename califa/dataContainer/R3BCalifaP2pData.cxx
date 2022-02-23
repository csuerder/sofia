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

#include "R3BCalifaP2pData.h"

R3BCalifaP2pData::R3BCalifaP2pData()
    : fEnergy1(NAN)
    , fEnergy2(NAN)
    , fPhi1(NAN)
    , fPhi2(NAN)
    , fTheta1(NAN)
    , fTheta2(NAN)
    , fTime(0)
    , fThetaSum(0.)
    , fPhiSum(0.)
    , fEnergySum(0.)
{
}

R3BCalifaP2pData::R3BCalifaP2pData(Double_t ene1, 
                                    Double_t ene2, 
                                    Double_t eneSum, 
                                    Double_t phi1, 
                                    Double_t phi2, 
                                    Double_t phiSum, 
                                    Double_t theta1, 
                                    Double_t theta2, 
                                    Double_t thetaSum, 
                                    ULong64_t time)
    : fEnergy1(ene1)
    , fEnergy2(ene2)
    , fPhi1(phi1)
    , fPhi2(phi2)
    , fTheta1(theta1)
    , fTheta2(theta2)
    , fTime(time)
    , fThetaSum(thetaSum)
    , fPhiSum(phiSum)
    , fEnergySum(eneSum)
{
}


ClassImp(R3BCalifaP2pData);
