#include "ReconstructionDataFormats/GlobalFwdTrack.h"
#include "DataFormatsMFT/TrackMFT.h"

#include<iostream>
#include<memory>

#include "ReconstructionDataFormats/TrackFwd.h"
#include "CommonConstants/MathConstants.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/Propagator.h"
#include "Field/MagneticField.h"
#include "SimulationDataFormat/TrackReference.h"
#include "SimulationDataFormat/MCEventHeader.h"
#include "SimulationDataFormat/MCTrack.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "SimulationDataFormat/BaseHits.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "TTree.h"
#include <TGeoGlobalMagField.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TEfficiency.h>

void MCTrackChecker_primary( const std::string o2sim_KineFile = "o2sim_Kine.root",
                     const std::string sig_KineFile = "sgn_Kine.root"){

  // Set MC track information
  std::unique_ptr<TFile> o2sim_KineFileIn(new TFile(o2sim_KineFile.c_str()));
  std::unique_ptr<TTree> o2SimKineTree((TTree*)o2sim_KineFileIn->Get("o2sim"));
  vector<o2::MCTrackT<float>>mcTrVec, *mcTr = &mcTrVec;
  o2SimKineTree->SetBranchAddress("MCTrack", &mcTr);
  o2::dataformats::MCEventHeader* eventHeader = nullptr;
  o2SimKineTree->SetBranchAddress("MCEventHeader.", &eventHeader);
  vector<o2::TrackReference>* trackRefs = nullptr;
  o2SimKineTree->SetBranchAddress("TrackRefs", &trackRefs);
  o2SimKineTree->GetEntry(0);

  //Parepare output objects
  TFile* output = new TFile("MCTrackinfo.root","recreate");
  std::unique_ptr<TH3F> MCTrackVertex(new TH3F("MCTrackVertex","",100,-0.2,0.2,100,-0.2,0.2,100,-0.2,0.2));
  std::unique_ptr<TH3F> MCTrackPosition(new TH3F("MCTrackPosition","",100,-0.2,0.2,100,-0.2,0.2,100,-0.2,0.2));
  std::unique_ptr<TH1F> MCTrackpT(new TH1F("MCTrackpT","",100,0,10));

  float vz_MC,P_MC;
  int pdgcode_MC,MotherId;
  int nEvent = o2SimKineTree->GetEntries();
  for (int iEvent=0; iEvent < nEvent; iEvent++) {

        //Get tree information of the "iTrack"th line
        o2SimKineTree->GetEntry(iEvent);
        int nMCTracks = (*mcTr).size();
        for (Int_t iMCTrack=0; iMCTrack<(*mcTr).size();++iMCTrack) {
                o2::MCTrackT<float> *thisTrack = &(*mcTr).at(iMCTrack);
                if(thisTrack->isPrimary()){

                        //Get Track information of the "iMCTrack"th line.
                        /*
                        vz_MC = thisTrack->GetStartVertexCoordinatesZ();//get particle production z-posiiton
                        P_MC = thisTrack->GetP();//get particle total momentum
                        pdgcode_MC = thisTrack->GetPdgCode();//get particle pdg code
                        MotherId = thisTrack->getMotherTrackId();
                        std::cout<<"Get Track information of the /iMCTrack/th line."<<std::endl;
                        std::cout<<"This Event = "<<iEvent<<std::endl;
                        std::cout<<"This Track = "<<iMCTrack<<std::endl;
                        std::cout<<"Vertex Z = "<<vz_MC<<std::endl;
                        std::cout<<"PDG code = "<<pdgcode_MC<<std::endl;
                        std::cout<<"Mother Id = "<<MotherId<<std::endl;

                        MotherId = thisTrack->getMotherTrackId();
                        //if(MotherId ==-1) continue;
                        while(MotherId !=-1){
                                //if(MotherId < iMCTrack) break;
                                thisTrack = &(mcTr->at(MotherId));
                                MotherId = thisTrack->getMotherTrackId();
                        }
                        //if(MotherId < iMCTrack) continue;
                        thisTrack = &(mcTr->at(iMCTrack));
                        */
                        vz_MC = thisTrack->GetStartVertexCoordinatesZ();//get particle production z-posiiton
                        P_MC = thisTrack->GetP();//get particle total momentum
                        pdgcode_MC = thisTrack->GetPdgCode();//get particle pdg code
                        MotherId = thisTrack->getMotherTrackId();
                        std::cout<<"Primary information"<<std::endl;
                        std::cout<<"This Event = "<<iEvent<<std::endl;
                        std::cout<<"This Track = "<<iMCTrack<<std::endl;
                        std::cout<<"Vertex Z = "<<vz_MC<<std::endl;
                        std::cout<<"PDG code = "<<pdgcode_MC<<std::endl;
                        std::cout<<"Mother Id = "<<MotherId<<std::endl;
                }
                     //MC track information
                auto vx_MC = thisTrack->GetStartVertexCoordinatesX();//get particle production x-posiiton
                auto vy_MC = thisTrack->GetStartVertexCoordinatesY();//get particle production y-posiiton
                vz_MC = thisTrack->GetStartVertexCoordinatesZ();//get particle production z-posiiton
                auto Pt_MC = thisTrack->GetPt();//get particle pt
                P_MC = thisTrack->GetP();//get particle total momentum
                pdgcode_MC = thisTrack->GetPdgCode();//get particle pdg code

                //Select only true muon track. 13 is muon pdg code
                // if(fabs(pdgcode_MC) != 211) {
                //   iTrack++;
                //   continue;
                // }

                //Fill data to object
                MCTrackVertex->Fill(vx_MC,vy_MC,vz_MC);
                MCTrackpT->Fill(Pt_MC);
                 }
        iEvent++;
         }
  //Save histgram
  output->WriteTObject(MCTrackVertex.get());
  output->WriteTObject(MCTrackPosition.get());
  output->WriteTObject(MCTrackpT.get());
}
