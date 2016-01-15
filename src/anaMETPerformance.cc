#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>
#include <stdio.h>

#include "UserCode/diall/interface/anaMETPerformance.h"
#include "UserCode/diall/interface/diParticle.h"
#include "UserCode/diall/interface/genParticle.h"
#include "UserCode/diall/interface/lwJetContainer.h"
#include "UserCode/diall/interface/particleBase.h"
#include "UserCode/diall/interface/pfParticle.h"


#include "TLorentzVector.h"

#include "TClass.h"

ClassImp(anaMETPerformance)

anaMETPerformance::anaMETPerformance(const char *name, const char *title) 
:anaBaseTask(name,title),
  fCheckPid(kFALSE),
  fMetType(),      
  fMinPt(0.),      
  fMuonsName(""),  
  fParticlesName(),
  fMuons(0x0),     
  fParticles(0x0), 
  fZs(0x0),
  fZsName(""),
  fMETPerformanceInfo(0x0),
  fNZCands(0),
  fMassZCands(0.),
  fPtZCands(0.),
  fEtaZCands(0.),
  fMET(0.),
  fuParaZll(0.),
  fuParaZllPt(0.),
  fuPerpZll(0.)
  
{
  
}

//----------------------------------------------------------
void anaMETPerformance::Exec(Option_t * /*option*/)
{

   fNZCands = 0;
  
   TLorentzVector met ; TLorentzVector l;
   int count = 0; 
   
   //if(!SelectEvent()) return;

   if(!fInitOutput) CreateOutputObjects();

   //Get objects from event
   if(!fMuons && !fMuonsName.IsNull()) {
     fMuons = dynamic_cast<TClonesArray*>(fEventObjects->FindObject(fMuonsName.Data()));
   }
   if(!fMuons) return;
   //Make array for Z candidates
   if(!fEventObjects->FindObject(fZsName) && !fZs) {
      fZs = new TClonesArray("diParticle");
      fZs->SetName(fZsName);
      fEventObjects->Add(fZs);
    }
   if(fZs) fZs->Delete();
   
   //Double_t cent = fHiEvent->GetCentrality();
   Int_t nmuons = fMuons->GetEntriesFast();
   Printf("nmuons: %d",nmuons);
   //   fh1NMuons->Fill(nmuons);
   if(nmuons<2) return;
   printf("anaMETPerformance executing\n"); 


   for (int i = 0; i < fMuons->GetEntriesFast(); i++) {
     particleBase *mu1 = static_cast<particleBase*>(fMuons->At(i));
     if(!mu1) {
       Printf("%s ERROR: couldn't get muon",GetName());
       continue;
     }
     if(fCheckPid)
       if(!CheckPid(mu1)) continue;
     count=0;
     for (int j = i+1; j < fMuons->GetEntriesFast(); j++) {
       particleBase *mu2 = static_cast<particleBase*>(fMuons->At(j));
       if(!mu2) {
         Printf("%s ERROR: couldn't get muon",GetName());
         continue;
       }
       TLorentzVector l1 = mu1->GetLorentzVector();
       TLorentzVector l2 = mu2->GetLorentzVector();
       TLorentzVector dimu = l1 + l2;
       
       //muons should be of opposite sign
       if(mu1->GetCharge()*mu2->GetCharge()<0) {

         if(fCheckPid)
           if(!CheckPid(mu2)) continue;
                   
         //fh3CentPtInvMass->Fill(cent,dimu.Pt(),dimu.M());
         
         //Store Z candidates in event
         if(fZs) {
	   diParticle *pPart = new ((*fZs)[fNZCands])
             diParticle(dimu.Pt(),
                        dimu.Eta(),
                        dimu.Phi(),
                        dimu.M(),
                        11);
           pPart->SetCharge(0);
           pPart->AddParticle(mu1);
           pPart->AddParticle(mu2);
	   count++;fNZCands++;
	   Printf("muon loop prin %d\n", fMuons->GetEntriesFast()-1);
         }
       } else {
         //fh3CentPtInvMassSC->Fill(cent,dimu.Pt(),dimu.M());
	 //Printf("muon loop prin %d\n", fMuons->GetEntriesFast()-1);
	 //flag_charge = false;
       }
             
     }//muon 2 loop
   }//muon 1 loop
   /*
   //Get particles from which MET will be calculated
   if(!fParticles && !fParticlesName.IsNull()) {
     //fEventObjects->Print();
     fParticles = dynamic_cast<TClonesArray*>(fEventObjects->FindObject(fParticlesName.Data()));
     if(!fParticles) {
       //check if in jet branch
       lwJetContainer *jetsCont = dynamic_cast<lwJetContainer*>(fEventObjects->FindObject(fParticlesName.Data()));
       if(jetsCont) fParticles = jetsCont->GetJets();
     }
   }
   
   if(!fParticles) {
     Printf("%s: WARNING: Couldn't locate %s branch",GetName(),fParticlesName.Data());
     return;
   }
   */
   //Get particles from which MET will be calculated
   
   if(!fParticles && !fParticlesName.IsNull()) {
     fParticles = dynamic_cast<TClonesArray*>(fEventObjects->FindObject(fParticlesName.Data()));
   }
   //Printf("%s: ",fParticlesName.Data());
   if(!fParticles) {
     Printf("%s: WARNING: Couldn't locate %s branch",GetName(),fParticlesName.Data());
     return;
   }
   
   //Get jet container
   //if(!fJetsCont && !fJetsName.IsNull())
   //  fJetsCont = dynamic_cast<lwJetContainer*>(fEventObjects->FindObject(fJetsName.Data()));


   //Double_t cent = 5.;//fHiEvent->GetCentrality();
   TLorentzVector p4(0.,0.,0.,0.); 
   Double_t sumEt = 0.;

   const Int_t nptmin = 10;
   Double_t ptarr[nptmin] {0.,1.,2.,3.,10.,20.,30.,40.,50.,60.};
   TLorentzVector r4[nptmin];
   for(Int_t j = 0; j<nptmin; ++j)
     r4[j].SetPtEtaPhiM(0.,0.,0.,0.);

   for (int i = 0; i < fParticles->GetEntriesFast(); i++) {
     particleBase *p = static_cast<particleBase*>(fParticles->At(i));
     if(!p) {
       Printf("%s ERROR: couldn't get particle",GetName());
       continue;
     }

     if(fMetType==kGen || fMetType==kPFRaw) {
     //  if(p->Pt() < fMinPt)
       //  continue;
       l = p->GetLorentzVector();
     }
     else if(fMetType==kVS) {
       pfParticle *pf = dynamic_cast<pfParticle*>(p);
       if(!pf) {
         Printf("%s ERROR: couldn't cast particle to pfParticle",GetName());
         return;
       }
       l.SetPtEtaPhiM(pf->PtVS(),pf->Eta(),pf->Phi(),pf->M());
     }
     else if(fMetType==kPuppi) {
       pfParticle *pf = dynamic_cast<pfParticle*>(p);
       if(!pf) {
         Printf("%s ERROR: couldn't cast particle to pfParticle",GetName());
         return;
       }
       if (pf->GetPuppiWeight()!=0) printf("weight %f\n", pf->GetPuppiWeight());
       l = pf->GetPuppiWeight()*p->GetLorentzVector();
     }
     
     for(Int_t j = 0; j<nptmin; ++j) {
       if(l.Pt()>ptarr[j])
         r4[j]+=l;
     }
    
     if(l.Pt() < fMinPt) continue;
     //fh3PtEtaPhi->Fill(l.Pt(),l.Eta(),l.Phi());
     p4+=l;
     sumEt+=l.Et();
 
   }//particle loop
   //Printf("muon loop %d\n", fParticles->GetEntriesFast());
	
   met = -p4;
   fMET = met.Pt();
   //printf("!! mpika edp %f %f\n", met.Pt(), p4.Pt());
   //fh2MetCent->Fill(cent,met.Pt());
   //fh2SumEtCent->Fill(cent,sumEt);

   Int_t nhists = TMath::Min(nptmin,10);
   for(Int_t j = 0; j<nhists; ++j) {
     TLorentzVector met2 = -r4[j];
     //fh2MetCentPtMin[j]->Fill(cent,met2.Pt());
   }
   std::cout<<fNZCands<<std::endl;
   int error;
   //fh1uParaZllPt->Fill(sumEt);
   for(int i = 0; i<fNZCands; ++i) { 
     diParticle *pPart = (diParticle*)fZs->At(i);
     //if (pPart->Pt()>20 && pPart->M()>60 && pPart->M()<120)
     //std::pair<double, double> u = compHadronicRecoilProjU(pPart,met, error, count);
     std::pair<double, double> u = compMETProjU(pPart,met.Px(), met.Py(), error, count);
     //if (pPart->Pt()>48 && pPart->Pt()<60)
     //fh1uParaZllPt->Fill(std::get<0>(u)+pPart->Pt());
     fMassZCands = pPart->M();
     fPtZCands = pPart->Pt();
     fEtaZCands = pPart->Eta();
     fuParaZll = std::get<0>(u);
     fuParaZllPt = std::get<0>(u)+pPart->Pt();
     fuPerpZll = std::get<1>(u);
     std::cout<< fuParaZllPt <<std::endl;;
     fMETPerformanceInfo->Fill();
   }


}

//----------------------------------------------------------                                                                         
bool anaMETPerformance::CheckPid(particleBase *p) {
  //check if generated particle is muon                                                                                              
  genParticle *gp = dynamic_cast<genParticle*>(p);
  if(!gp) return kFALSE;
  if(abs(gp->GetPID())==13) return kTRUE;
  return kFALSE;
}

std::pair<double, double> 
anaMETPerformance::compMETProjU(diParticle* zP4, double metPx, double metPy, int& errorFlag, int count)
{
  if ( zP4->Pt() == 0. ) {
    Warning ("compMEtProjU", " Failed to compute projection, because Z0 candidate has zero Pt --> returning dummy solution !!");
    errorFlag = 1;
    return std::pair<double, double>(0., 0.);
  }
  
  double qX = zP4->Px();
  double qY = zP4->Py();
  double qT = TMath::Sqrt(qX*qX + qY*qY);
  
  double uX = -metPx;
  double uY = -metPy;
  uX -= qX;
  uY -= qY;
  
  
  double u1 = (uX*qX + uY*qY)/qT;
  double u2 = (uX*qY - uY*qX)/qT;

  return std::pair<double, double>(u1,u2);
}

std::pair<double, double> 
anaMETPerformance::compHadronicRecoilProjU(diParticle* zP4, TLorentzVector MET, int& errorFlag, int count)
{
  if ( zP4->Pt() == 0. ) {
    Warning ("compMEtProjU", " Failed to compute projection, because Z0 candidate has zero Pt --> returning dummy solution !!");
    errorFlag = 1;
    return std::pair<double, double>(0., 0.);
  }
  
  double qX = zP4->Px();
  double qY = zP4->Py();
  double qT = TMath::Sqrt(qX*qX + qY*qY);
  
  TLorentzVector zBoson = zP4->GetLorentzVector();
  TLorentzVector hadronicRecoil = -(MET+zBoson);

  double uX = hadronicRecoil.Px();
  double uY = hadronicRecoil.Py();
  //uX -= qX;
  //uY -= qY;
  
  double u1 = (uX*qX + uY*qY)/qT;
  double u2 = (uX*qY - uY*qX)/qT;

  return std::pair<double, double>(u1,u2);
}

//----------------------------------------------------------
void 
anaMETPerformance::CreateOutputObjects() {

  anaBaseTask::CreateOutputObjects();

  if(!fOutput) {
    Printf("anaMETPerformance: fOutput not present");
    return;
  }

   fMETPerformanceInfo = new TTree("METPerformanceInf","Tree with info related to MET performance with recoils");
   fMETPerformanceInfo->Branch("nZCands", &fNZCands, "nZCands/I");
   fMETPerformanceInfo->Branch("MassZCands",&fMassZCands, "MassZCands/F");
   fMETPerformanceInfo->Branch("PtZCands",&fPtZCands, "PtZCands/F");
   fMETPerformanceInfo->Branch("EtaZCands",&fEtaZCands, "EtaZCands/F");
   fMETPerformanceInfo->Branch("MET",&fMET, "MET/F");
   fMETPerformanceInfo->Branch("uParaZll",&fuParaZll, "uParaZll/F");
   fMETPerformanceInfo->Branch("uParaZllPt",&fuParaZllPt, "uParaZllPt/F");
   fMETPerformanceInfo->Branch("uPerpZllPt",&fuPerpZll, "uPerpZll/F");
   
   fOutput->Add(fMETPerformanceInfo);

  //for(Int_t j = 0; j<10; ++j) {
  //  fh2MetCentPtMin[j] = new TH2F(Form("fh2MetCentPtMin%d",j),"fh2MetCent;centrality;MET",100,0,100,500,0,1000.);
  //  fOutput->Add(fh2MetCentPtMin[j]);
  // }
    
  
}

