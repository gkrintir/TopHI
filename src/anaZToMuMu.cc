#include "UserCode/diall/interface/anaZToMuMu.h"
#include "UserCode/diall/interface/genParticle.h"
#include "UserCode/diall/interface/diParticle.h"

#include "TLorentzVector.h"

#include "TClass.h"

ClassImp(anaZToMuMu)

//----------------------------------------------------------
anaZToMuMu::anaZToMuMu(const char *name, const char *title) 
:anaBaseTask(name,title),
  fMuonsName(""),
  fMuons(0x0),
  fZsName(""),
  fZs(0x0),
  fCheckPid(kFALSE),
  fh1NMuons(),
  fh3CentPtInvMass(),
  fh3CentPtInvMassSC()
{

}

//----------------------------------------------------------
void anaZToMuMu::Exec(Option_t * /*option*/)
{
   anaBaseTask::Exec();
   //printf("anaZToMuMu executing !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! \n");

   if(!SelectEvent()) return;
  
   if(!fInitOutput) CreateOutputObjects();

   //Get objects from event

   //Get particles from which MET will be calculated
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

   Double_t cent = fHiEvent->GetCentrality();
   Int_t nmuons = fMuons->GetEntriesFast();
   Printf("nmuons: %d",nmuons);
   fh1NMuons->Fill(nmuons);
   if(nmuons<2) return;

   Printf("muon loop %d\n", fMuons->GetEntriesFast()-1);
   for (int i = 0; i < fMuons->GetEntriesFast()-1; i++) {
     particleBase *mu1 = static_cast<particleBase*>(fMuons->At(i));
     if(!mu1) {
       Printf("%s ERROR: couldn't get muon",GetName());
       continue;
     }
     if(fCheckPid)
       if(!CheckPid(mu1)) continue;
     //Printf("i: %d",i);
     int count = 0;
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
       if(mu1->GetCharge()*mu2->GetCharge()>0) {

         if(fCheckPid)
           if(!CheckPid(mu2)) continue;
                   
         fh3CentPtInvMass->Fill(cent,dimu.Pt(),dimu.M());
         
         //Store Z candidates in event
         if(fZs) {
           diParticle *pPart = new ((*fZs)[count])
             diParticle(dimu.Pt(),
                        dimu.Eta(),
                        dimu.Phi(),
                        dimu.M(),
                        11);
           pPart->SetCharge(0);
           pPart->AddParticle(mu1);
           pPart->AddParticle(mu2);
           ++count;
         }
       } else {
         fh3CentPtInvMassSC->Fill(cent,dimu.Pt(),dimu.M());
       }
             
     }//muon 2 loop
   }//muon 1 loop

}

//----------------------------------------------------------
Bool_t anaZToMuMu::CheckPid(particleBase *p) {
  //check if generated particle is muon
  genParticle *gp = dynamic_cast<genParticle*>(p);
  if(!gp) return kFALSE;
  if(abs(gp->GetPID())==13) return kTRUE;
  return kFALSE;
}

//----------------------------------------------------------
void anaZToMuMu::CreateOutputObjects() {

  anaBaseTask::CreateOutputObjects();

  if(!fOutput) {
    Printf("anaZToMuMu: fOutput not present");
    return;
  }

  fh1NMuons = new TH1F("fh1NMuons","fh1NMuons;#it{N}_{muons}",100,0,100.);
  fOutput->Add(fh1NMuons);
  
  fh3CentPtInvMass = new TH3F("fh3CentPtInvMass","fh3CentPtInvMass;centrality;#it{p}_{T};#it{M}_{inv}",100,0,100,200,0,200.,200,0.,200.);
  fOutput->Add(fh3CentPtInvMass);

  fh3CentPtInvMassSC = new TH3F("fh3CentPtInvMassSC","fh3CentPtInvMassSC;centrality;#it{p}_{T};#it{M}_{inv}",100,0,100,200,0,200.,200,0.,200.);
  fOutput->Add(fh3CentPtInvMassSC);

}
