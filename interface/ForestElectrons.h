#ifndef ForestElectrons_h
#define ForestElectrons_h

#define maxForestElectrons 10

#include <iostream>
#include <vector>
#include "TBranch.h"

class ForestElectrons {
public :
   ForestElectrons(){};
   ~ForestElectrons(){};

   // Declaration of leaf types
   // Event info
   Int_t           Run;
   Int_t           Event;
   Int_t           Lumi;
   Int_t           CentBin;
   Float_t         vx;
   Float_t         vy;
   Float_t         vz;
   // GenParticle info
   Int_t           Gen_nptl;
   Int_t           Gen_pid[maxForestElectrons];   //[Gen_nptl]
   Int_t           Gen_mom[maxForestElectrons];   //[Gen_nptl] pid mother
   Int_t           Gen_status[maxForestElectrons];   //[Gen_nptl]
   Float_t         Gen_p[maxForestElectrons];   //[Gen_nptl]
   Float_t         Gen_pt[maxForestElectrons];   //[Gen_nptl]
   Float_t         Gen_eta[maxForestElectrons];   //[Gen_nptl]
   Float_t         Gen_phi[maxForestElectrons];   //[Gen_nptl]
   // RecoElectron info
   Int_t           nEle;
   std::vector<Int_t>     *eleCharge;
   std::vector<Int_t>     *eleChargeConsistent;
   std::vector<Float_t>   *eleEn;
   std::vector<Float_t>   *eleD0;
   std::vector<Float_t>   *eleDz;
   std::vector<Float_t>   *eleD0Err;
   std::vector<Float_t>   *eleDzErr;
   std::vector<Float_t>   *eleTrkPt;
   std::vector<Float_t>   *eleTrkEta;
   std::vector<Float_t>   *eleTrkPhi;
   std::vector<Int_t>     *eleTrkCharge;
   std::vector<Float_t>   *eleTrkChi2;
   std::vector<Float_t>   *eleTrkNdof;
   std::vector<Float_t>   *eleTrkNormalizedChi2;
   std::vector<Int_t>     *eleTrkValidHits;
   std::vector<Int_t>     *eleTrkLayers;
   std::vector<Float_t>   *elePt;
   std::vector<Float_t>   *eleEta;
   std::vector<Float_t>   *elePhi;
   std::vector<Float_t>   *eleSCEn;
   std::vector<Float_t>   *eleESEn;
   std::vector<Float_t>   *eleSCEta;
   std::vector<Float_t>   *eleSCPhi;
   std::vector<Float_t>   *eleSCRawEn;
   std::vector<Float_t>   *eleSCEtaWidth;
   std::vector<Float_t>   *eleSCPhiWidth;
   std::vector<Float_t>   *eleHoverE;
   std::vector<Float_t>   *eleEoverP;
   std::vector<Float_t>   *eleEoverPInv;
   std::vector<Float_t>   *eleBrem;
   std::vector<Float_t>   *eledEtaAtVtx;
   std::vector<Float_t>   *eledPhiAtVtx;
   std::vector<Float_t>   *eleSigmaIEtaIEta;
   std::vector<Float_t>   *eleSigmaIEtaIEta_2012;
   std::vector<Float_t>   *eleSigmaIPhiIPhi;
   std::vector<Int_t>     *eleMissHits;
   std::vector<Float_t>   *eleESEffSigmaRR;
   std::vector<Float_t>   *elePFChIso;
   std::vector<Float_t>   *elePFPhoIso;
   std::vector<Float_t>   *elePFNeuIso;
   std::vector<Float_t>   *elePFPUIso;
   std::vector<Float_t>   *elePFChIso03;
   std::vector<Float_t>   *elePFPhoIso03;
   std::vector<Float_t>   *elePFNeuIso03;
   std::vector<Float_t>   *elePFChIso04;
   std::vector<Float_t>   *elePFPhoIso04;
   std::vector<Float_t>   *elePFNeuIso04;
   std::vector<Float_t>   *eleBC1E;
   std::vector<Float_t>   *eleBC1Eta;
   std::vector<Float_t>   *eleBC2E;
   std::vector<Float_t>   *eleBC2Eta;
   // Dielectron
   Int_t                        Di_npair;
   std::vector<Float_t>         *Di_vProb=0;   //[Di_npair]
   std::vector<Float_t>         *Di_mass=0;   //[Di_npair]
   std::vector<Float_t>         *Di_e=0;   //[Di_npair]
   std::vector<Float_t>         *Di_pt=0;   //[Di_npair]
   std::vector<Float_t>         *Di_pt1=0;   //[Di_npair]
   std::vector<Float_t>         *Di_pt2=0;   //[Di_npair]
   std::vector<Float_t>         *Di_eta=0;   //[Di_npair]
   std::vector<Float_t>         *Di_eta1=0;   //[Di_npair]
   std::vector<Float_t>         *Di_eta2=0;   //[Di_npair]
   std::vector<Float_t>         *Di_phi=0;   //[Di_npair]
   std::vector<Float_t>         *Di_phi1=0;   //[Di_npair]
   std::vector<Float_t>         *Di_phi2=0;   //[Di_npair]
   std::vector<Int_t>           *Di_charge1=0;   //[Di_npair]
   std::vector<Int_t>           *Di_charge2=0;   //[Di_npair]
   std::vector<Int_t>           *Di_isArb1=0;   //[Di_npair]
   std::vector<Int_t>           *Di_isArb2=0;   //[Di_npair]
   std::vector<Float_t>         *Di_nTrkHit1=0;   //[Di_npair]
   std::vector<Float_t>         *Di_nTrkHit2=0;   //[Di_npair]
   std::vector<Float_t>         *Di_trkChi2_1=0;   //[Di_npair]
   std::vector<Float_t>         *Di_trkChi2_2=0;   //[Di_npair]
   std::vector<Float_t>         *Di_glbChi2_1=0;   //[Di_npair]
   std::vector<Float_t>         *Di_glbChi2_2=0;   //[Di_npair]
   std::vector<Float_t>         *Di_dxy1=0;   //[Di_npair]
   std::vector<Float_t>         *Di_dxy2=0;   //[Di_npair]
   std::vector<Float_t>         *Di_dz1=0;   //[Di_npair]
   std::vector<Float_t>         *Di_dz2=0;   //[Di_npair]

   // List of branches
   // Event info
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_cbin;   //!
   TBranch        *b_vx;   //!
   TBranch        *b_vy;   //!
   TBranch        *b_vz;   //!
   // GenParticle info
   TBranch        *b_Gen_nptl;   //!
   TBranch        *b_Gen_pid;   //!
   TBranch        *b_Gen_mom;   //!
   TBranch        *b_Gen_status;   //!
   TBranch        *b_Gen_p;   //!
   TBranch        *b_Gen_pt;   //!
   TBranch        *b_Gen_eta;   //!
   TBranch        *b_Gen_phi;   //!
   // RecoElectron info
   TBranch        *b_nEle;
   TBranch        *b_eleCharge;
   TBranch        *b_eleChargeConsistent;
   TBranch        *b_eleEn;
   TBranch        *b_eleD0;
   TBranch        *b_eleDz;
   TBranch        *b_eleD0Err;
   TBranch        *b_eleDzErr;
   TBranch        *b_eleTrkPt;
   TBranch        *b_eleTrkEta;
   TBranch        *b_eleTrkPhi;
   TBranch        *b_eleTrkCharge;
   TBranch        *b_eleTrkChi2;
   TBranch        *b_eleTrkNdof;
   TBranch        *b_eleTrkNormalizedChi2;
   TBranch        *b_eleTrkValidHits;
   TBranch        *b_eleTrkLayers;
   TBranch        *b_elePt;
   TBranch        *b_eleEta;
   TBranch        *b_elePhi;
   TBranch        *b_eleSCEn;
   TBranch        *b_eleESEn;
   TBranch        *b_eleSCEta;
   TBranch        *b_eleSCPhi;
   TBranch        *b_eleSCRawEn;
   TBranch        *b_eleSCEtaWidth;
   TBranch        *b_eleSCPhiWidth;
   TBranch        *b_eleHoverE;
   TBranch        *b_eleEoverP;
   TBranch        *b_eleEoverPInv;
   TBranch        *b_eleBrem;
   TBranch        *b_eledEtaAtVtx;
   TBranch        *b_eledPhiAtVtx;
   TBranch        *b_eleSigmaIEtaIEta;
   TBranch        *b_eleSigmaIEtaIEta_2012;
   TBranch        *b_eleSigmaIPhiIPhi;
   TBranch        *b_eleMissHits;
   TBranch        *b_eleESEffSigmaRR;
   TBranch        *b_elePFChIso;
   TBranch        *b_elePFPhoIso;
   TBranch        *b_elePFNeuIso;
   TBranch        *b_elePFPUIso;
   TBranch        *b_elePFChIso03;
   TBranch        *b_elePFPhoIso03;
   TBranch        *b_elePFNeuIso03;
   TBranch        *b_elePFChIso04;
   TBranch        *b_elePFPhoIso04;
   TBranch        *b_elePFNeuIso04;
   TBranch        *b_eleBC1E;
   TBranch        *b_eleBC1Eta;
   TBranch        *b_eleBC2E;
   TBranch        *b_eleBC2Eta;
   //Dielectrons
   TBranch        *b_Di_npair;   //!
   TBranch        *b_Di_vProb;   //!
   TBranch        *b_Di_mass;   //!
   TBranch        *b_Di_e;   //!
   TBranch        *b_Di_pt;   //!
   TBranch        *b_Di_pt1;   //!
   TBranch        *b_Di_pt2;   //!
   TBranch        *b_Di_eta;   //!
   TBranch        *b_Di_eta1;   //!
   TBranch        *b_Di_eta2;   //!
   TBranch        *b_Di_phi;   //!
   TBranch        *b_Di_phi1;   //!
   TBranch        *b_Di_phi2;   //!
   TBranch        *b_Di_charge1;   //!
   TBranch        *b_Di_charge2;   //!
   TBranch        *b_Di_isArb1;   //!
   TBranch        *b_Di_isArb2;   //!
   TBranch        *b_Di_nTrkHit1;   //!
   TBranch        *b_Di_nTrkHit2;   //!
   TBranch        *b_Di_trkChi2_1;   //!
   TBranch        *b_Di_trkChi2_2;   //!
   TBranch        *b_Di_glbChi2_1;   //!
   TBranch        *b_Di_glbChi2_2;   //!
   TBranch        *b_Di_dxy1;   //!
   TBranch        *b_Di_dxy2;   //!
   TBranch        *b_Di_dz1;   //!
   TBranch        *b_Di_dz2;   //!
};
#endif
