//
// producer of particle flow candidates
//

#include "UserCode/diall/interface/lwJetFromForestProducer.h"
#include "UserCode/diall/interface/lwJet.h"

ClassImp(lwJetFromForestProducer)

//__________________________________________________________
lwJetFromForestProducer::lwJetFromForestProducer() :
inputBase("lwJetFromForestProducer"),
  flwJetContName("flwJetCont"),
  flwJetContainer(0x0),
  flwGenJetContName("flwGenJetCont"),
  flwGenJetContainer(0x0),
  fForestJets(),
  fRadius(-1.)
{
  //default constructor
}

//__________________________________________________________
lwJetFromForestProducer::lwJetFromForestProducer(const char *name) :
  inputBase(name),
  flwJetContName("flwJetCont"),
  flwJetContainer(0x0),
  flwGenJetContName("flwGenJetCont"),
  flwGenJetContainer(0x0),
  fForestJets(),
  fRadius(-1.)
{
  //standard constructor
}

//__________________________________________________________
Bool_t lwJetFromForestProducer::Init() {

  if(!inputBase::Init()) return kFALSE;
  
  if(fInputMode==hiForest) {
    if (fChain->GetBranch("nref"))
      fChain->SetBranchAddress("nref", &fForestJets.nref, &fForestJets.b_nref);
    if (fChain->GetBranch("rawpt"))
      fChain->SetBranchAddress("rawpt", fForestJets.rawpt, &fForestJets.b_rawpt);
    if (fChain->GetBranch("jtpt"))
      fChain->SetBranchAddress("jtpt", fForestJets.jtpt, &fForestJets.b_jtpt);
    if (fChain->GetBranch("jteta"))
      fChain->SetBranchAddress("jteta", fForestJets.jteta, &fForestJets.b_jteta);
    if (fChain->GetBranch("jty"))
      fChain->SetBranchAddress("jty", fForestJets.jty, &fForestJets.b_jty);
    if (fChain->GetBranch("jtphi"))
      fChain->SetBranchAddress("jtphi", fForestJets.jtphi, &fForestJets.b_jtphi);
    if (fChain->GetBranch("jtpu"))
      fChain->SetBranchAddress("jtpu", fForestJets.jtpu, &fForestJets.b_jtpu);
    if (fChain->GetBranch("jtm"))
      fChain->SetBranchAddress("jtm", fForestJets.jtm, &fForestJets.b_jtm);
    if (fChain->GetBranch("discr_ssvHighEff"))
      fChain->SetBranchAddress("discr_ssvHighEff", fForestJets.discr_ssvHighEff, &fForestJets.b_discr_ssvHighEff);
    if (fChain->GetBranch("discr_ssvHighPur"))
      fChain->SetBranchAddress("discr_ssvHighPur", fForestJets.discr_ssvHighPur, &fForestJets.b_discr_ssvHighPur);
    if (fChain->GetBranch("discr_csvMva"))
      fChain->SetBranchAddress("discr_csvMva", fForestJets.discr_csvMva, &fForestJets.b_discr_csvMva);
    if (fChain->GetBranch("discr_csvSimple"))
      fChain->SetBranchAddress("discr_csvSimple", fForestJets.discr_csvSimple, &fForestJets.b_discr_csvSimple);
    if (fChain->GetBranch("refparton_flavor"))
      fChain->SetBranchAddress("refparton_flavor", fForestJets.refparton_flavor, &fForestJets.b_refparton_flavor);
    if (fChain->GetBranch("refparton_flavorForB"))
      fChain->SetBranchAddress("refparton_flavorForB", fForestJets.refparton_flavorForB, &fForestJets.b_refparton_flavorForB);

    if (fChain->GetBranch("ngen"))
      fChain->SetBranchAddress("ngen", &fForestJets.ngen, &fForestJets.b_ngen);
    if (fChain->GetBranch("genpt"))
      fChain->SetBranchAddress("genpt", fForestJets.genpt, &fForestJets.b_genpt);
    if (fChain->GetBranch("geneta"))
      fChain->SetBranchAddress("geneta", fForestJets.geneta, &fForestJets.b_geneta);
    if (fChain->GetBranch("geny"))
      fChain->SetBranchAddress("geny", fForestJets.geny, &fForestJets.b_geny);
    if (fChain->GetBranch("genphi"))
      fChain->SetBranchAddress("genphi", fForestJets.genphi, &fForestJets.b_genphi);
    if (fChain->GetBranch("gensubid"))
      fChain->SetBranchAddress("gensubid", fForestJets.gensubid, &fForestJets.b_gensubid);

    fInit = kTRUE;
  }
  return kTRUE;
}

//__________________________________________________________
Bool_t lwJetFromForestProducer::InitEventObjects() {
  //Create event objects
  if(!fEventObjects) {
    Printf("%s: fEventObjects does not exist. Cannot store output",GetName());
    return kFALSE;
  } else {
    if(!fEventObjects->FindObject(flwJetContName) && !flwJetContName.IsNull()) {
      flwJetContainer = new lwJetContainer(flwJetContName);
      flwJetContainer->Init();
      flwJetContainer->SetJetRadius(fRadius);
      fEventObjects->Add(flwJetContainer);
    }
    if(!fEventObjects->FindObject(flwGenJetContName) && !flwGenJetContName.IsNull()) {
      flwGenJetContainer = new lwJetContainer(flwGenJetContName);
      flwGenJetContainer->Init();
      flwGenJetContainer->SetJetRadius(fRadius);
      fEventObjects->Add(flwGenJetContainer);
    }
  }

  return kTRUE;
}

//__________________________________________________________
Bool_t lwJetFromForestProducer::Run(Long64_t entry) {

  //run analysis
  Long64_t centry = LoadTree(entry);
  if(centry<0) return kFALSE;

  if(!InitEventObjects()) return kFALSE;
 
  //clear array
  flwJetContainer->ClearVec();
  
  //put jets of this event in array
  Int_t jetCount = 0;
  //Printf("%s: njets: %d",GetName(),fForestJets.nref);
  for(Int_t i = 0; i<fForestJets.nref; i++) {
    lwJet *jet = new lwJet(fForestJets.jtpt[i],
                           fForestJets.jteta[i],
                           fForestJets.jtphi[i],
                           fForestJets.jtm[i]);
    jet->SetRefToParton(fForestJets.refparton_flavor[i]);
    jet->SetRefToPartonForB(fForestJets.refparton_flavorForB[i]);
    jet->SetCsvSimpleDiscr(fForestJets.discr_csvSimple[i]);
    flwJetContainer->AddJet(jet,jetCount);
    ++jetCount;
  }
  //Printf("jetCount: %d",jetCount);
  flwJetContainer->SortJets();

  //only produce gen jets if requested and if they exist
  if(flwGenJetContainer && fChain->GetBranch("ngen")) {
    Printf("doing gen jets");
    if(fForestJets.ngen>0) {
      flwGenJetContainer->ClearVec();
      Int_t genJetCount = 0;
      for(Int_t i = 0; i<fForestJets.ngen; i++) {
        lwJet *genjet = new lwJet(fForestJets.genpt[i],
                                  fForestJets.geneta[i],
                                  fForestJets.genphi[i],
                                  0.);
        flwGenJetContainer->AddJet(genjet,genJetCount);
        ++genJetCount;
      }
      flwGenJetContainer->SortJets();
    }
  }
  
  return kTRUE;
}

//__________________________________________________________
Long64_t lwJetFromForestProducer::LoadTree(Long64_t entry) {

  //overloaded LoadTree function 
  if(!fChain) {
    Printf("fChain doesn't exist");
    return -5;
  }
  
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
    Init();
    //  Printf("%lld fCurrent: %d",entry,fCurrent);
  }

  // fChain->SetMakeClass(1);
 
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) {
    Printf("lwJetFromForestProducer: centry smaller than 0");
    return centry;  
  } 

  fChain->GetEntry(entry);
  
  return centry;
}
